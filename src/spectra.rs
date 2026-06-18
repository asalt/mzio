use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::fs::OpenOptions;
use std::io::Stdout;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{mpsc as std_mpsc, Mutex, OnceLock};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use crossterm::{
    event::{
        self, DisableMouseCapture, EnableMouseCapture, Event as CrosstermEvent, KeyCode, KeyEvent,
    },
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use mzdata::{
    io::{DetailLevel, MzMLReader, SpectrumSource},
    prelude::*,
    spectrum::SignalContinuity,
};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    widgets::{
        canvas::{Canvas, Line as CanvasLine},
        Block, Borders, List, ListItem, Paragraph, Wrap,
    },
    Frame, Terminal,
};
use tokio::{
    sync::mpsc,
    time::{self, Duration},
};

mod backend;
use backend::{IndexResult, MzmlFileInfo, SpectrumData, SpectrumMeta, SpectrumOwned};

type CrosstermTerminal = Terminal<CrosstermBackend<Stdout>>;

const MAX_CACHE: usize = 64;
const MAX_CACHE_BYTES: usize = 256 * 1024 * 1024;
const MAX_DRAW_PEAKS: usize = 800;
const MAX_PLOT_BINS: usize = 50_000;
const SVG_WIDTH: u32 = 1200;
const SVG_HEIGHT: u32 = 600;
const SVG_BINS: usize = 4000;

#[derive(Clone, Debug)]
pub struct IndexCacheOptions {
    pub enabled: bool,
    pub refresh: bool,
    pub path: Option<PathBuf>,
}

impl Default for IndexCacheOptions {
    fn default() -> Self {
        Self {
            enabled: true,
            refresh: false,
            path: None,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct SpectraOptions {
    pub index_cache: IndexCacheOptions,
}

static SPECTRA_DEBUG_IDX: OnceLock<u32> = OnceLock::new();
static SPECTRA_DEBUG_FILE: OnceLock<Option<Mutex<fs::File>>> = OnceLock::new();

fn spectra_debug_enabled() -> bool {
    std::env::var_os("MZIO_SPECTRA_DEBUG").is_some()
}

fn spectra_debug_target(idx: u32) -> bool {
    if !spectra_debug_enabled() {
        return false;
    }
    if let Ok(val) = std::env::var("MZIO_SPECTRA_DEBUG_IDX") {
        if let Ok(target) = val.parse::<u32>() {
            return target == idx;
        }
    }
    *SPECTRA_DEBUG_IDX.get_or_init(|| idx) == idx
}

fn spectra_debug_log(msg: impl AsRef<str>) {
    if !spectra_debug_enabled() {
        return;
    }

    let file = SPECTRA_DEBUG_FILE.get_or_init(|| {
        let _ = fs::create_dir_all("logs");
        match OpenOptions::new()
            .create(true)
            .append(true)
            .open("logs/spectra_debug.log")
        {
            Ok(file) => Some(Mutex::new(file)),
            Err(_) => None,
        }
    });

    let Some(file) = file.as_ref() else {
        return;
    };

    let ts_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_millis();
    let mut file = file.lock().expect("debug log mutex poisoned");
    let _ = writeln!(file, "[{ts_ms}] {}", msg.as_ref());
}

#[derive(Debug)]
enum SpectraEvent {
    Tick,
    Key(KeyEvent),
    IndexLoaded {
        metas: Vec<SpectrumMeta>,
        file_info: MzmlFileInfo,
        cached: bool,
        cache_path: Option<PathBuf>,
        cache_warning: Option<String>,
    },
    DataLoaded {
        idx: u32,
        data: SpectrumOwned,
        primary: bool,
    },
    LoadFailed {
        idx: u32,
        primary: bool,
        error: String,
    },
    Status(String),
}

#[derive(Debug)]
enum LoaderCmd {
    SetOffsets { offsets: Vec<u64> },
    Load { idx: u32, primary: bool },
    Quit,
}

#[derive(Debug, Clone, Copy)]
struct LoadingState {
    idx: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PlotMode {
    Auto,
    Sticks,
    Line,
}

impl PlotMode {
    fn next(self) -> Self {
        match self {
            PlotMode::Auto => PlotMode::Sticks,
            PlotMode::Sticks => PlotMode::Line,
            PlotMode::Line => PlotMode::Auto,
        }
    }

    fn label(self) -> &'static str {
        match self {
            PlotMode::Auto => "auto",
            PlotMode::Sticks => "sticks",
            PlotMode::Line => "line",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PlotRenderMode {
    Sticks,
    Line,
}

impl PlotRenderMode {
    fn label(self) -> &'static str {
        match self {
            PlotRenderMode::Sticks => "sticks",
            PlotRenderMode::Line => "line",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Zoom {
    Full,
    Low,
    Mid,
    High,
    Precursor,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RightPane {
    Spectrum,
    Chromatogram,
    Map,
    Dashboard,
}

impl RightPane {
    fn toggle(self) -> Self {
        match self {
            RightPane::Spectrum => RightPane::Chromatogram,
            RightPane::Chromatogram => RightPane::Map,
            RightPane::Map => RightPane::Dashboard,
            RightPane::Dashboard => RightPane::Spectrum,
        }
    }

    fn label(self) -> &'static str {
        match self {
            RightPane::Spectrum => "spectrum",
            RightPane::Chromatogram => "chrom",
            RightPane::Map => "map",
            RightPane::Dashboard => "overview",
        }
    }
}

impl Zoom {
    fn next(self) -> Self {
        match self {
            Zoom::Full => Zoom::Low,
            Zoom::Low => Zoom::Mid,
            Zoom::Mid => Zoom::High,
            Zoom::High => Zoom::Precursor,
            Zoom::Precursor => Zoom::Full,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Zoom::Full => "full",
            Zoom::Low => "0-500",
            Zoom::Mid => "500-1000",
            Zoom::High => "1000-2000",
            Zoom::Precursor => "precursor ±75",
        }
    }
}

struct SpectraApp {
    path: PathBuf,
    metas: Vec<SpectrumMeta>,
    file_info: Option<MzmlFileInfo>,
    selected: usize,
    list_offset: usize,
    cache: HashMap<u32, SpectrumOwned>,
    cache_order: VecDeque<u32>,
    cache_bytes: usize,
    normalize: bool,
    zoom: Zoom,
    status: String,
    loading: Option<LoadingState>,
    pending_primary: Option<u32>,
    loader_tx: std_mpsc::Sender<LoaderCmd>,
    pending_count: Option<usize>,
    search_mode: bool,
    search_query: String,
    last_search: Option<String>,
    last_search_forward: bool,
    plot_mode: PlotMode,
    plot_cache: Option<PlotCache>,
    chrom_cache: Option<ChromCache>,
    map_cache: Option<MapCache>,
    right_pane: RightPane,
}

impl SpectraApp {
    fn new(path: PathBuf, loader_tx: std_mpsc::Sender<LoaderCmd>) -> Self {
        Self {
            path,
            metas: Vec::new(),
            file_info: None,
            selected: 0,
            list_offset: 0,
            cache: HashMap::new(),
            cache_order: VecDeque::new(),
            cache_bytes: 0,
            normalize: true,
            zoom: Zoom::Full,
            status: String::from("Loading mzML metadata..."),
            loading: None,
            pending_primary: None,
            loader_tx,
            pending_count: None,
            search_mode: false,
            search_query: String::new(),
            last_search: None,
            last_search_forward: true,
            plot_mode: PlotMode::Auto,
            plot_cache: None,
            chrom_cache: None,
            map_cache: None,
            right_pane: RightPane::Spectrum,
        }
    }

    fn set_index(
        &mut self,
        metas: Vec<SpectrumMeta>,
        file_info: MzmlFileInfo,
        cached: bool,
        cache_path: Option<PathBuf>,
        cache_warning: Option<String>,
    ) {
        let mut loader_ok = true;
        if !metas.is_empty() {
            let offsets = metas.iter().map(|m| m.offset).collect::<Vec<_>>();
            loader_ok = self
                .loader_tx
                .send(LoaderCmd::SetOffsets { offsets })
                .is_ok();
        }
        self.metas = metas;
        self.file_info = Some(file_info);
        self.selected = 0;
        self.list_offset = 0;
        let mut notes = Vec::<String>::new();
        if let Some(warn) = cache_warning.as_ref() {
            if !warn.trim().is_empty() {
                notes.push(warn.trim().to_string());
            }
        } else if !cached && cache_path.is_some() {
            notes.push("cache saved".to_string());
        }
        if !loader_ok {
            notes.push("loader stopped".to_string());
        }

        let mut status = if cached {
            format!("Loaded cached index for {} spectra", self.metas.len())
        } else {
            format!("Indexed {} spectra", self.metas.len())
        };
        if !notes.is_empty() {
            status.push_str(&format!(" ({})", notes.join("; ")));
        }
        self.status = status;
        self.plot_cache = None;
        self.chrom_cache = None;
        self.map_cache = None;
    }

    fn set_status(&mut self, msg: impl Into<String>) {
        self.status = msg.into();
    }

    fn toggle_normalize(&mut self) {
        self.normalize = !self.normalize;
        self.plot_cache = None;
        self.chrom_cache = None;
        self.map_cache = None;
    }

    fn cycle_zoom(&mut self) {
        self.zoom = self.zoom.next();
        self.plot_cache = None;
    }

    fn cycle_plot_mode(&mut self) {
        self.plot_mode = self.plot_mode.next();
        self.plot_cache = None;
    }

    fn toggle_right_pane(&mut self) {
        self.right_pane = self.right_pane.toggle();
    }

    fn move_selection(&mut self, delta: isize) {
        if self.metas.is_empty() {
            return;
        }
        let len = self.metas.len() as isize;
        let mut new_idx = self.selected as isize + delta;
        if new_idx < 0 {
            new_idx = 0;
        } else if new_idx >= len {
            new_idx = len - 1;
        }
        self.selected = new_idx as usize;
    }

    fn jump_to(&mut self, idx: usize) {
        if idx < self.metas.len() {
            self.selected = idx;
        }
    }

    fn ensure_visible(&mut self, viewport: usize) {
        if self.selected < self.list_offset {
            self.list_offset = self.selected;
        } else if self.selected >= self.list_offset.saturating_add(viewport) {
            self.list_offset = self.selected.saturating_sub(viewport.saturating_sub(1));
        }
    }

    fn current_meta(&self) -> Option<&SpectrumMeta> {
        self.metas.get(self.selected)
    }

    fn has_data(&self, idx: u32) -> bool {
        self.cache.contains_key(&idx)
    }

    fn data(&self, idx: u32) -> Option<SpectrumData<'_>> {
        self.cache.get(&idx).map(|owned| SpectrumData {
            mz: Cow::Borrowed(&*owned.mz),
            intensity: Cow::Borrowed(&*owned.intensity),
        })
    }

    fn insert_data(&mut self, idx: u32, data: SpectrumOwned) {
        if self.cache.contains_key(&idx) {
            return;
        }
        let bytes = data.approx_bytes();
        self.cache.insert(idx, data);
        self.cache_order.push_back(idx);
        self.cache_bytes = self.cache_bytes.saturating_add(bytes);
        while self.cache_order.len() > MAX_CACHE
            || (self.cache_bytes > MAX_CACHE_BYTES && self.cache_order.len() > 1)
        {
            let Some(evicted) = self.cache_order.pop_front() else {
                break;
            };
            if let Some(evicted_data) = self.cache.remove(&evicted) {
                self.cache_bytes = self.cache_bytes.saturating_sub(evicted_data.approx_bytes());
            }
            if matches!(self.plot_cache, Some(ref cache) if cache.idx == evicted) {
                self.plot_cache = None;
            }
        }
        if matches!(self.plot_cache, Some(ref cache) if cache.idx == idx) {
            self.plot_cache = None;
        }
    }

    fn maybe_start_pending_primary(&mut self) -> bool {
        if self.loading.is_some() {
            return false;
        }
        let Some(idx) = self.pending_primary.take() else {
            return false;
        };
        if self.has_data(idx) {
            return false;
        }
        self.request_load(idx, true);
        true
    }

    fn request_load(&mut self, idx: u32, primary: bool) {
        if self.has_data(idx) {
            if primary {
                self.set_status(format!("Cached spectrum {idx}"));
            }
            return;
        }
        if let Some(loading) = self.loading {
            if loading.idx == idx {
                return;
            }
            if primary {
                self.pending_primary = Some(idx);
                self.set_status(format!(
                    "Queued spectrum {idx} (loading {}...)",
                    loading.idx
                ));
            }
            return;
        }
        self.loading = Some(LoadingState { idx });
        if primary {
            self.set_status(format!("Loading spectrum {idx}..."));
        }
        if self
            .loader_tx
            .send(LoaderCmd::Load { idx, primary })
            .is_err()
        {
            self.loading = None;
            self.set_status("mzML loader thread stopped".to_string());
        }
    }

    fn clear_count(&mut self) {
        self.pending_count = None;
    }

    fn take_count(&mut self) -> usize {
        let count = self.pending_count.take().unwrap_or(1);
        count.max(1)
    }

    fn push_count_digit(&mut self, ch: char) {
        if let Some(d) = ch.to_digit(10) {
            let current = self.pending_count.unwrap_or(0);
            let next = current.saturating_mul(10).saturating_add(d as usize).max(1);
            self.pending_count = Some(next);
        }
    }

    fn start_search(&mut self) {
        self.search_mode = true;
        self.search_query.clear();
        self.clear_count();
        self.set_status("Search: ");
    }

    fn cancel_search(&mut self) {
        self.search_mode = false;
        self.search_query.clear();
    }

    fn apply_search(&mut self, forward: bool) -> bool {
        if self.metas.is_empty() {
            return false;
        }

        let pattern = if !self.search_query.is_empty() {
            self.search_query.clone()
        } else if let Some(ref last) = self.last_search {
            last.clone()
        } else {
            return false;
        };

        let needle = pattern.to_ascii_lowercase();
        self.last_search = Some(pattern.clone());
        self.last_search_forward = forward;

        let len = self.metas.len();
        if len == 0 {
            return false;
        }

        let mut idx = self.selected;
        for _ in 0..len {
            idx = if forward {
                (idx + 1) % len
            } else {
                (idx + len - 1) % len
            };

            let Some(meta) = self.metas.get(idx).cloned() else {
                continue;
            };

            let haystack = format!(
                "{} ms{} rt={:?}min {} prec={:?} z={:?} points={}",
                meta.scan_id,
                meta.ms_level,
                meta.rt_minutes,
                meta.continuity,
                meta.precursor_mz,
                meta.charge,
                meta.points
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            )
            .to_ascii_lowercase();

            if haystack.contains(&needle) {
                self.selected = idx;
                self.ensure_visible(10);
                self.set_status(format!("Search '{}' → idx {}", pattern, idx));
                self.request_load(meta.idx, true);
                return true;
            }
        }

        self.set_status(format!("Search '{}' not found", pattern));
        false
    }
}

#[derive(Debug, Clone)]
struct PlotCache {
    idx: u32,
    zoom: Zoom,
    normalize: bool,
    mode: PlotRenderMode,
    canvas_width: u16,
    x_bounds: (f64, f64),
    points: Vec<(f64, f64)>,
    y_max: f64,
}

#[derive(Debug, Clone)]
struct ChromCache {
    normalize: bool,
    canvas_width: u16,
    rt_bounds: (f64, f64),
    tic_points: Vec<(f64, f64)>,
    tic_y_max: f64,
    bpc_points: Vec<(f64, f64)>,
    bpc_y_max: f64,
}

#[derive(Debug, Clone)]
struct MapCache {
    normalize: bool,
    canvas_size: (u16, u16),
    rt_bounds: (f64, f64),
    mz_bounds: (f64, f64),
    max_intensity: f64,
    buckets: [Vec<(f64, f64)>; 5],
}

pub async fn run_spectra_demo(path: PathBuf, options: SpectraOptions) -> anyhow::Result<()> {
    if spectra_debug_enabled() {
        spectra_debug_log(format!("spectra demo start: {}", path.display()));
    }
    if let Ok(val) = std::env::var("MZIO_SPECTRA_DEBUG_DUMP_IDX") {
        let idx = val
            .parse::<u32>()
            .with_context(|| format!("invalid MZIO_SPECTRA_DEBUG_DUMP_IDX={val:?}"))?;
        spectra_debug_log(format!("spectra debug dump: idx={idx}"));
        let mut reader = MzMLReader::open_path(path.as_path())
            .with_context(|| format!("failed to open mzML at {}", path.display()))?;
        reader.detail_level = DetailLevel::Lazy;
        let data = backend::load_spectrum(&mut reader, idx)?;
        spectra_debug_log(format!(
            "spectra debug dump complete: idx={idx} mz.len={} intensity.len={}",
            data.mz.len(),
            data.intensity.len()
        ));
        return Ok(());
    }
    let mut terminal = init_terminal()?;
    let result = run_spectra_loop(&mut terminal, path, options).await;
    restore_terminal(&mut terminal)?;
    result
}

async fn run_spectra_loop(
    terminal: &mut CrosstermTerminal,
    path: PathBuf,
    options: SpectraOptions,
) -> anyhow::Result<()> {
    let (tx, mut rx) = mpsc::channel::<SpectraEvent>(256);
    spawn_input_listener(tx.clone());
    spawn_tick(tx.clone(), Duration::from_millis(150));
    spawn_indexer(path.clone(), tx.clone(), options.index_cache.clone());
    let loader_tx = spawn_loader(path.clone(), tx.clone());

    let mut app = SpectraApp::new(path.clone(), loader_tx);
    terminal.draw(|f| draw(f, &mut app))?;

    while let Some(ev) = rx.recv().await {
        match ev {
            SpectraEvent::Tick => {
                terminal.draw(|f| draw(f, &mut app))?;
            }
            SpectraEvent::Key(key) => {
                if matches!(key.code, KeyCode::Char('q') | KeyCode::Esc) {
                    break;
                }
                handle_key(&mut app, key);
                terminal.draw(|f| draw(f, &mut app))?;
            }
            SpectraEvent::IndexLoaded {
                metas,
                file_info,
                cached,
                cache_path,
                cache_warning,
            } => {
                if spectra_debug_enabled() {
                    spectra_debug_log(format!(
                        "IndexLoaded: {} spectra (cached={cached})",
                        metas.len()
                    ));
                }
                app.set_index(metas, file_info, cached, cache_path, cache_warning);
                if let Some(meta) = app.current_meta() {
                    app.request_load(meta.idx, true);
                }
                terminal.draw(|f| draw(f, &mut app))?;
            }
            SpectraEvent::DataLoaded { idx, data, primary } => {
                if matches!(app.loading, Some(loading) if loading.idx == idx) {
                    app.loading = None;
                }
                let stats = data.stats;
                if spectra_debug_target(idx) {
                    spectra_debug_log(format!(
                        "DataLoaded idx={idx} primary={primary} points={} mz[min,max]=[{:.4},{:.4}] base_peak_mz={:.4} base_peak_intensity={:.4e}",
                        stats.points,
                        stats.mz_min,
                        stats.mz_max,
                        stats.base_peak_mz,
                        stats.base_peak_intensity
                    ));
                }
                app.insert_data(idx, data);
                if let Some(meta) = app.metas.get_mut(idx as usize) {
                    meta.points = Some(stats.points);
                    meta.mz_min = Some(stats.mz_min);
                    meta.mz_max = Some(stats.mz_max);
                    meta.base_peak_mz = Some(stats.base_peak_mz);
                    meta.base_peak_intensity = Some(stats.base_peak_intensity);
                }
                if primary {
                    app.set_status(format!("Loaded spectrum {idx}"));
                }
                let started_pending = app.maybe_start_pending_primary();
                if primary && !started_pending {
                    prefetch_neighbors(&mut app);
                }
                terminal.draw(|f| draw(f, &mut app))?;
            }
            SpectraEvent::LoadFailed {
                idx,
                primary,
                error,
            } => {
                if matches!(app.loading, Some(loading) if loading.idx == idx) {
                    app.loading = None;
                }
                if spectra_debug_target(idx) {
                    spectra_debug_log(format!("LoadFailed idx={idx} primary={primary}: {error}"));
                }
                if primary {
                    app.set_status(format!("load error: {error}"));
                } else {
                    app.set_status(format!("prefetch error: {error}"));
                }
                let _ = app.maybe_start_pending_primary();
                terminal.draw(|f| draw(f, &mut app))?;
            }
            SpectraEvent::Status(msg) => {
                app.set_status(msg);
                terminal.draw(|f| draw(f, &mut app))?;
            }
        }
    }

    let _ = app.loader_tx.send(LoaderCmd::Quit);

    Ok(())
}

fn handle_key(app: &mut SpectraApp, key: KeyEvent) {
    let mut moved = false;

    // When in search mode, interpret keys as editing the search query.
    if app.search_mode {
        match key.code {
            KeyCode::Esc => {
                app.cancel_search();
                app.set_status("Search cancelled");
            }
            KeyCode::Enter => {
                app.search_mode = false;
                let _ = app.apply_search(true);
            }
            KeyCode::Backspace => {
                app.search_query.pop();
                app.set_status(format!("Search: {}", app.search_query));
            }
            KeyCode::Char(c) => {
                app.search_query.push(c);
                app.set_status(format!("Search: {}", app.search_query));
            }
            _ => {}
        }
        return;
    }

    match key.code {
        KeyCode::Char(c) if c.is_ascii_digit() => {
            app.push_count_digit(c);
        }
        KeyCode::Char('/') => {
            app.start_search();
        }
        KeyCode::Char('s') => {
            app.clear_count();
            moved = app.apply_search(true);
        }
        KeyCode::Char('S') => {
            app.clear_count();
            moved = app.apply_search(false);
        }
        KeyCode::Up | KeyCode::Char('k') => {
            let count = app.take_count();
            app.move_selection(-(count as isize));
            moved = true;
        }
        KeyCode::Down | KeyCode::Char('j') => {
            let count = app.take_count();
            app.move_selection(count as isize);
            moved = true;
        }
        KeyCode::PageUp => {
            let count = app.take_count();
            app.move_selection(-20 * count as isize);
            moved = true;
        }
        KeyCode::PageDown => {
            let count = app.take_count();
            app.move_selection(20 * count as isize);
            moved = true;
        }
        KeyCode::Home | KeyCode::Char('g') => {
            app.clear_count();
            app.jump_to(0);
            moved = true;
        }
        KeyCode::End | KeyCode::Char('G') => {
            app.clear_count();
            if !app.metas.is_empty() {
                app.jump_to(app.metas.len().saturating_sub(1));
                moved = true;
            }
        }
        KeyCode::Char('n') => {
            app.clear_count();
            app.toggle_normalize();
        }
        KeyCode::Char('z') | KeyCode::Tab => {
            app.clear_count();
            app.cycle_zoom();
        }
        KeyCode::Char('m') => {
            app.clear_count();
            app.cycle_plot_mode();
        }
        KeyCode::Char('o') => {
            app.clear_count();
            app.toggle_right_pane();
        }
        KeyCode::Char('p') => {
            app.clear_count();
            match export_current_spectrum_svg(app) {
                Ok(path) => app.set_status(format!("Exported SVG: {}", path.display())),
                Err(err) => app.set_status(format!("Export failed: {err:?}")),
            }
        }
        KeyCode::Char('P') => {
            app.clear_count();
            match export_current_spectrum_png(app) {
                Ok(path) => app.set_status(format!("Exported PNG: {}", path.display())),
                Err(err) => app.set_status(format!("Export failed: {err:?}")),
            }
        }
        KeyCode::Enter => {
            app.clear_count();
            if let Some(meta) = app.current_meta() {
                app.request_load(meta.idx, true);
            }
        }
        _ => {
            app.clear_count();
        }
    }

    if moved {
        prefetch_neighbors(app);
    }
}

fn spawn_tick(tx: mpsc::Sender<SpectraEvent>, interval: Duration) {
    tokio::spawn(async move {
        let mut ticker = time::interval(interval);
        loop {
            ticker.tick().await;
            if tx.send(SpectraEvent::Tick).await.is_err() {
                break;
            }
        }
    });
}

fn spawn_input_listener(tx: mpsc::Sender<SpectraEvent>) {
    tokio::spawn(async move {
        loop {
            match event::poll(Duration::from_millis(50)) {
                Ok(true) => match event::read() {
                    Ok(CrosstermEvent::Key(key)) => {
                        if tx.send(SpectraEvent::Key(key)).await.is_err() {
                            break;
                        }
                    }
                    Ok(_) => {}
                    Err(_) => break,
                },
                Ok(false) => {}
                Err(_) => break,
            }
        }
    });
}

fn spawn_indexer(path: PathBuf, tx: mpsc::Sender<SpectraEvent>, cache: IndexCacheOptions) {
    tokio::spawn(async move {
        let progress_tx = tx.clone();
        let result = tokio::task::spawn_blocking(move || {
            backend::index_mzml(path.as_path(), &cache, |msg| {
                let _ = progress_tx.blocking_send(SpectraEvent::Status(msg));
            })
        })
        .await;
        match result {
            Ok(Ok(IndexResult {
                metas,
                file_info,
                cached,
                cache_path,
                cache_warning,
            })) => {
                let _ = tx
                    .send(SpectraEvent::IndexLoaded {
                        metas,
                        file_info,
                        cached,
                        cache_path,
                        cache_warning,
                    })
                    .await;
            }
            Ok(Err(err)) => {
                let _ = tx
                    .send(SpectraEvent::Status(format!("index error: {err:?}")))
                    .await;
            }
            Err(err) => {
                let _ = tx
                    .send(SpectraEvent::Status(format!("task error: {err:?}")))
                    .await;
            }
        }
    });
}

fn spawn_loader(path: PathBuf, tx: mpsc::Sender<SpectraEvent>) -> std_mpsc::Sender<LoaderCmd> {
    let (cmd_tx, cmd_rx) = std_mpsc::channel::<LoaderCmd>();
    std::thread::spawn(move || {
        let mut reader = match fs::File::open(path.clone()) {
            Ok(file) => match MzMLReader::open_file(file) {
                Ok(mut reader) => {
                    reader.detail_level = DetailLevel::Lazy;
                    Some(reader)
                }
                Err(err) => {
                    let _ = tx.blocking_send(SpectraEvent::Status(format!(
                        "failed to initialize mzML reader: {err:?}"
                    )));
                    None
                }
            },
            Err(err) => {
                let _ = tx.blocking_send(SpectraEvent::Status(format!(
                    "failed to open mzML for loading: {err:?}"
                )));
                None
            }
        };
        let mut offsets_ready = false;

        while let Ok(cmd) = cmd_rx.recv() {
            match cmd {
                LoaderCmd::Quit => break,
                LoaderCmd::SetOffsets { offsets } => {
                    let Some(reader) = reader.as_mut() else {
                        continue;
                    };
                    let mut index = mzdata::io::OffsetIndex::new("spectrum".to_string());
                    for (i, offset) in offsets.iter().enumerate() {
                        index.insert(i.to_string(), *offset);
                    }
                    index.init = true;
                    reader.set_index(index);
                    offsets_ready = true;
                }
                LoaderCmd::Load { idx, primary } => {
                    let Some(reader) = reader.as_mut() else {
                        let _ = tx.blocking_send(SpectraEvent::LoadFailed {
                            idx,
                            primary,
                            error: "mzML loader not initialized".to_string(),
                        });
                        continue;
                    };

                    if !offsets_ready {
                        let _ = tx.blocking_send(SpectraEvent::LoadFailed {
                            idx,
                            primary,
                            error: "mzML loader offset index not initialized".to_string(),
                        });
                        continue;
                    }

                    let res = backend::load_spectrum(reader, idx);
                    match res {
                        Ok(data) => {
                            let _ =
                                tx.blocking_send(SpectraEvent::DataLoaded { idx, data, primary });
                        }
                        Err(err) => {
                            let _ = tx.blocking_send(SpectraEvent::LoadFailed {
                                idx,
                                primary,
                                error: format!("{err:?}"),
                            });
                        }
                    }
                }
            }
        }
    });
    cmd_tx
}

fn prefetch_neighbors(app: &mut SpectraApp) {
    if app.loading.is_some() || app.pending_primary.is_some() || app.metas.is_empty() {
        return;
    }
    let len = app.metas.len();
    let mut candidates = Vec::new();
    let span = 4usize;
    for delta in 1..=span {
        if let Some(idx) = app.selected.checked_add(delta) {
            if idx < len {
                candidates.push(idx);
            }
        }
        if delta <= app.selected {
            candidates.push(app.selected - delta);
        }
    }

    for idx in candidates {
        if let Some(meta) = app.metas.get(idx) {
            if !app.has_data(meta.idx) {
                app.request_load(meta.idx, false);
                break;
            }
        }
    }
}

fn draw(frame: &mut Frame<'_>, app: &mut SpectraApp) {
    let size = frame.size();
    let columns = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Length(42), Constraint::Min(20)])
        .split(size);

    draw_list(frame, app, columns[0]);
    draw_main(frame, app, columns[1]);
}

fn draw_list(frame: &mut Frame<'_>, app: &mut SpectraApp, area: Rect) {
    let mut items: Vec<ListItem> = Vec::new();
    let view_height = area.height.saturating_sub(2) as usize;
    let total = app.metas.len();

    app.ensure_visible(view_height.max(1));
    let offset = app.list_offset.min(total.saturating_sub(1));

    for (i, meta) in app
        .metas
        .iter()
        .enumerate()
        .skip(offset)
        .take(view_height.max(1))
    {
        let selected = i == app.selected;
        let rt_str = meta
            .rt_minutes
            .map(|v| format!("{v:6.2}m"))
            .unwrap_or_else(|| "   -  ".to_string());
        let prec = meta
            .precursor_mz
            .map(|v| format!("{v:7.2}"))
            .unwrap_or_else(|| "  -   ".to_string());
        let ch = meta
            .charge
            .map(|c| format!("{c:+}"))
            .unwrap_or_else(|| "-".to_string());
        let sig = continuity_short(meta.continuity);
        let points = meta
            .points
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string());
        let line = format!(
            "{:>5} | ms{} | {} | rt {} | prec {} | z {} | pts {} | {}",
            meta.idx, meta.ms_level, sig, rt_str, prec, ch, points, meta.scan_id
        );
        let mut item = ListItem::new(line);
        if selected {
            item = item.style(Style::default().fg(Color::Yellow));
        }
        items.push(item);
    }

    let list = List::new(items)
        .block(Block::default().title("Spectra").borders(Borders::ALL))
        .highlight_style(Style::default().fg(Color::Yellow))
        .highlight_symbol("➜ ");

    frame.render_widget(list, area);
}

fn draw_main(frame: &mut Frame<'_>, app: &mut SpectraApp, area: Rect) {
    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(7),
            Constraint::Min(10),
            Constraint::Length(1),
        ])
        .split(area);

    draw_meta_panel(frame, app, rows[0]);
    match app.right_pane {
        RightPane::Spectrum => draw_plot(frame, app, rows[1]),
        RightPane::Chromatogram => draw_chromatogram(frame, app, rows[1]),
        RightPane::Map => draw_map(frame, app, rows[1]),
        RightPane::Dashboard => draw_dashboard(frame, app, rows[1]),
    }

    let status = Paragraph::new(app.status.as_str())
        .style(Style::default().fg(Color::DarkGray))
        .block(Block::default().borders(Borders::TOP));
    frame.render_widget(status, rows[2]);
}

fn draw_meta_panel(frame: &mut Frame<'_>, app: &SpectraApp, area: Rect) {
    let mut lines = Vec::new();
    lines.push(format!(
        "File: {}",
        app.path
            .file_name()
            .map(|s| s.to_string_lossy())
            .unwrap_or_else(|| app.path.to_string_lossy())
    ));
    if let Some(meta) = app.current_meta() {
        lines.push(format!(
            "Scan {} | ms{} | rt={}min | {} | precursor={} m/z | z={} | pts={}",
            meta.scan_id,
            meta.ms_level,
            meta.rt_minutes
                .map(|v| format!("{v:.2}"))
                .unwrap_or_else(|| "-".to_string()),
            meta.continuity,
            meta.precursor_mz
                .map(|v| format!("{v:.4}"))
                .unwrap_or_else(|| "-".to_string()),
            meta.charge
                .map(|v| format!("{v}"))
                .unwrap_or_else(|| "-".to_string()),
            meta.points
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string())
        ));
        if let (Some(bp_mz), Some(bp_i)) = (meta.base_peak_mz, meta.base_peak_intensity) {
            lines.push(format!("Base peak: {:.4} m/z @ {:.3e}", bp_mz, bp_i));
        }
    } else {
        lines.push("No spectra loaded yet.".to_string());
    }

    lines.push(format!(
        "[Enter] load | o: view={} | n: norm={} | z/Tab: zoom={} | m: mode={} | p: SVG | P: PNG | /,s,S: search | q: quit{}",
        app.right_pane.label(),
        if app.normalize { "on" } else { "off" },
        app.zoom.label(),
        app.plot_mode.label(),
        app.pending_count
            .map(|c| format!("  (count {})", c))
            .unwrap_or_default()
    ));

    if app.search_mode {
        lines.push(format!("Search: /{}", app.search_query));
    } else if let Some(ref q) = app.last_search {
        let dir = if app.last_search_forward {
            "↓"
        } else {
            "↑"
        };
        lines.push(format!("Search: /{} ({})", q, dir));
    }

    let widget = Paragraph::new(lines.join("\n"))
        .wrap(Wrap { trim: false })
        .block(Block::default().borders(Borders::ALL).title("Info"));
    frame.render_widget(widget, area);
}

fn draw_dashboard(frame: &mut Frame<'_>, app: &SpectraApp, area: Rect) {
    let mut lines: Vec<String> = Vec::new();

    let total = app.metas.len();
    let mut ms1 = 0usize;
    let mut ms2_plus = 0usize;
    let mut rt_min = f32::INFINITY;
    let mut rt_max = -f32::INFINITY;
    let mut has_rt = false;
    let mut prec_count = 0usize;
    let mut prec_min = f64::INFINITY;
    let mut prec_max = -f64::INFINITY;

    for meta in app.metas.iter() {
        if meta.ms_level <= 1 {
            ms1 = ms1.saturating_add(1);
        } else {
            ms2_plus = ms2_plus.saturating_add(1);
        }
        if let Some(rt) = meta.rt_minutes {
            has_rt = true;
            rt_min = rt_min.min(rt);
            rt_max = rt_max.max(rt);
        }
        if let Some(pmz) = meta.precursor_mz {
            prec_count = prec_count.saturating_add(1);
            prec_min = prec_min.min(pmz);
            prec_max = prec_max.max(pmz);
        }
    }

    lines.push(format!(
        "Indexed: {total} spectra  (MS1 {ms1}, MS2+ {ms2_plus})"
    ));
    if has_rt {
        lines.push(format!("RT range: {rt_min:.2}–{rt_max:.2} min"));
    }
    if prec_count > 0 && prec_min.is_finite() && prec_max.is_finite() {
        lines.push(format!(
            "Precursors: {prec_count}  (m/z {prec_min:.2}–{prec_max:.2})"
        ));
    }

    if let Some(info) = app.file_info.as_ref() {
        if info.run_id.is_some() || info.start_time.is_some() {
            lines.push(String::new());
            let run_id = info
                .run_id
                .as_deref()
                .filter(|s| !s.trim().is_empty())
                .unwrap_or("-");
            let start = info
                .start_time
                .as_deref()
                .filter(|s| !s.trim().is_empty())
                .unwrap_or("-");
            lines.push(format!("Run id: {run_id}"));
            lines.push(format!("Start time: {start}"));
        }

        if info.default_instrument_id.is_some() || info.spectrum_count_hint.is_some() {
            let inst = info
                .default_instrument_id
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string());
            let hint = info
                .spectrum_count_hint
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string());
            lines.push(format!("Spectrum count hint: {hint}"));
            lines.push(format!("Default instrument id: {inst}"));
        }

        if !info.contents.is_empty() {
            lines.push(String::new());
            lines.push(format!("Contents: {}", info.contents.join(", ")));
        }

        if !info.source_files.is_empty() {
            lines.push(String::new());
            lines.push("Source files:".to_string());
            for item in info.source_files.iter() {
                lines.push(format!("  {item}"));
            }
        }

        if !info.instrument_summaries.is_empty() {
            lines.push(String::new());
            lines.push("Instrument configurations:".to_string());
            for item in info.instrument_summaries.iter() {
                lines.push(format!("  {item}"));
            }
        }

        if !info.software_summaries.is_empty() {
            lines.push(String::new());
            lines.push("Software:".to_string());
            for item in info.software_summaries.iter() {
                lines.push(format!("  {item}"));
            }
        }
    } else {
        lines.push(String::new());
        lines.push("Metadata: loading...".to_string());
    }

    let widget = Paragraph::new(lines.join("\n"))
        .wrap(Wrap { trim: false })
        .block(Block::default().borders(Borders::ALL).title("Overview"));
    frame.render_widget(widget, area);
}

fn draw_chromatogram(frame: &mut Frame<'_>, app: &mut SpectraApp, area: Rect) {
    if app.metas.is_empty() {
        let msg = Paragraph::new("Indexing spectra...")
            .block(Block::default().borders(Borders::ALL).title("Chromatogram"));
        frame.render_widget(msg, area);
        return;
    }

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Min(4),
            Constraint::Min(4),
            Constraint::Length(1),
        ])
        .split(area);
    let tic_area = rows[0];
    let bpc_area = rows[1];
    let axis_area = rows[2];

    let current_rt = app.current_meta().and_then(|m| m.rt_minutes).map(f64::from);
    let Some(cache) = chromatogram_cached(app, tic_area.width.max(1)) else {
        let msg = Paragraph::new("No retention time metadata available")
            .block(Block::default().borders(Borders::ALL).title("Chromatogram"));
        frame.render_widget(msg, area);
        return;
    };

    if cache.tic_points.is_empty() {
        let msg = Paragraph::new("No TIC metadata (MS:1000285) found")
            .block(Block::default().borders(Borders::ALL).title("TIC (MS1)"));
        frame.render_widget(msg, tic_area);
    } else {
        let y_max = cache.tic_y_max.max(1e-6);
        let canvas = Canvas::default()
            .block(Block::default().borders(Borders::ALL).title("TIC (MS1)"))
            .x_bounds([cache.rt_bounds.0, cache.rt_bounds.1])
            .y_bounds([0.0, y_max])
            .paint(|ctx| {
                for w in cache.tic_points.windows(2) {
                    let (x1, y1) = w[0];
                    let (x2, y2) = w[1];
                    ctx.draw(&CanvasLine {
                        x1,
                        y1,
                        x2,
                        y2,
                        color: Color::LightGreen,
                    });
                }
                if let Some(rt) = current_rt {
                    if rt >= cache.rt_bounds.0 && rt <= cache.rt_bounds.1 {
                        ctx.draw(&CanvasLine {
                            x1: rt,
                            y1: 0.0,
                            x2: rt,
                            y2: y_max,
                            color: Color::DarkGray,
                        });
                    }
                }
            });
        frame.render_widget(canvas, tic_area);
    }

    if cache.bpc_points.is_empty() {
        let msg = Paragraph::new("No base peak intensity metadata (MS:1000505) found")
            .block(Block::default().borders(Borders::ALL).title("BPC (MS1)"));
        frame.render_widget(msg, bpc_area);
    } else {
        let y_max = cache.bpc_y_max.max(1e-6);
        let canvas = Canvas::default()
            .block(Block::default().borders(Borders::ALL).title("BPC (MS1)"))
            .x_bounds([cache.rt_bounds.0, cache.rt_bounds.1])
            .y_bounds([0.0, y_max])
            .paint(|ctx| {
                for w in cache.bpc_points.windows(2) {
                    let (x1, y1) = w[0];
                    let (x2, y2) = w[1];
                    ctx.draw(&CanvasLine {
                        x1,
                        y1,
                        x2,
                        y2,
                        color: Color::LightCyan,
                    });
                }
                if let Some(rt) = current_rt {
                    if rt >= cache.rt_bounds.0 && rt <= cache.rt_bounds.1 {
                        ctx.draw(&CanvasLine {
                            x1: rt,
                            y1: 0.0,
                            x2: rt,
                            y2: y_max,
                            color: Color::DarkGray,
                        });
                    }
                }
            });
        frame.render_widget(canvas, bpc_area);
    }

    let mid = (cache.rt_bounds.0 + cache.rt_bounds.1) / 2.0;
    let axis_label = format!(
        "{:.2}    {:.2}    {:.2} min",
        cache.rt_bounds.0, mid, cache.rt_bounds.1
    );
    let axis = Paragraph::new(axis_label)
        .style(Style::default().fg(Color::DarkGray))
        .alignment(Alignment::Center);
    frame.render_widget(axis, axis_area);
}

fn draw_map(frame: &mut Frame<'_>, app: &mut SpectraApp, area: Rect) {
    use ratatui::symbols::Marker;
    use ratatui::widgets::canvas::Points;

    if app.metas.is_empty() {
        let msg = Paragraph::new("Indexing spectra...")
            .block(Block::default().borders(Borders::ALL).title("Map"));
        frame.render_widget(msg, area);
        return;
    }

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Min(5), Constraint::Length(1)])
        .split(area);
    let canvas_area = rows[0];
    let axis_area = rows[1];

    let normalize = app.normalize;
    let current = app.current_meta().and_then(|m| {
        if m.ms_level > 1 {
            return None;
        }
        Some((m.rt_minutes? as f64, m.base_peak_mz?))
    });

    let Some(cache) = base_peak_map_cached(app, canvas_area) else {
        let msg = Paragraph::new("No MS1 base peak metadata available").block(
            Block::default()
                .borders(Borders::ALL)
                .title("MS1 Base Peak Map"),
        );
        frame.render_widget(msg, area);
        return;
    };

    const COLORS: [Color; 5] = [
        Color::DarkGray,
        Color::Blue,
        Color::Cyan,
        Color::Yellow,
        Color::Red,
    ];

    let title = format!(
        "MS1 Base Peak Map (max {:.3e})",
        cache.max_intensity.max(0.0)
    );
    let canvas = Canvas::default()
        .block(Block::default().borders(Borders::ALL).title(title))
        .marker(Marker::Braille)
        .x_bounds([cache.rt_bounds.0, cache.rt_bounds.1])
        .y_bounds([cache.mz_bounds.0, cache.mz_bounds.1])
        .paint(|ctx| {
            for (coords, color) in cache.buckets.iter().zip(COLORS.iter()) {
                if !coords.is_empty() {
                    ctx.draw(&Points {
                        coords,
                        color: *color,
                    });
                }
            }
            if let Some((rt, mz)) = current {
                let point = [(rt, mz)];
                ctx.draw(&Points {
                    coords: &point,
                    color: Color::White,
                });
            }
        });
    frame.render_widget(canvas, canvas_area);

    let x_mid = (cache.rt_bounds.0 + cache.rt_bounds.1) / 2.0;
    let axis_label = format!(
        "{:.2}    {:.2}    {:.2} min   |   m/z {:.0}–{:.0}   |   n: {}",
        cache.rt_bounds.0,
        x_mid,
        cache.rt_bounds.1,
        cache.mz_bounds.0,
        cache.mz_bounds.1,
        if normalize {
            "log colors"
        } else {
            "linear colors"
        }
    );
    let axis = Paragraph::new(axis_label)
        .style(Style::default().fg(Color::DarkGray))
        .alignment(Alignment::Center);
    frame.render_widget(axis, axis_area);
}

fn draw_plot(frame: &mut Frame<'_>, app: &mut SpectraApp, area: Rect) {
    let Some(meta) = app.current_meta() else {
        let msg = Paragraph::new("Indexing spectra...")
            .block(Block::default().borders(Borders::ALL).title("Spectrum"));
        frame.render_widget(msg, area);
        return;
    };

    let owned = app.cache.get(&meta.idx).cloned();
    if owned.is_none() {
        let mut msg = "Press Enter to load this spectrum".to_string();
        if let Some(loading) = app.loading {
            if loading.idx == meta.idx {
                msg = format!("Loading spectrum {} ...", meta.scan_id);
            }
        }
        let widget =
            Paragraph::new(msg).block(Block::default().borders(Borders::ALL).title("Spectrum"));
        frame.render_widget(widget, area);
        return;
    }
    let owned = owned.unwrap();
    let data = SpectrumData {
        mz: Cow::Borrowed(&*owned.mz),
        intensity: Cow::Borrowed(&*owned.intensity),
    };

    let x_bounds = compute_bounds(app.zoom, meta, &data);

    // Split plot area into the main canvas and a one-line x-axis label row.
    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Min(3), Constraint::Length(1)])
        .split(area);
    let canvas_area = rows[0];
    let axis_area = rows[1];

    let render_mode =
        resolve_plot_mode(app.plot_mode, meta.continuity, owned.stats.points as usize);
    let (points, y_max) =
        plot_points_cached(app, meta.idx, &data, x_bounds, render_mode, canvas_area);

    let canvas = Canvas::default()
        .block(Block::default().borders(Borders::ALL).title("Spectrum"))
        .x_bounds([x_bounds.0, x_bounds.1])
        .y_bounds([0.0, y_max.max(1.0)])
        .paint(|ctx| match render_mode {
            PlotRenderMode::Sticks => {
                for (x, y) in points.iter() {
                    ctx.draw(&CanvasLine {
                        x1: *x,
                        y1: 0.0,
                        x2: *x,
                        y2: *y,
                        color: Color::LightCyan,
                    });
                }
            }
            PlotRenderMode::Line => {
                for w in points.windows(2) {
                    let (x1, y1) = w[0];
                    let (x2, y2) = w[1];
                    ctx.draw(&CanvasLine {
                        x1,
                        y1,
                        x2,
                        y2,
                        color: Color::LightCyan,
                    });
                }
            }
        });

    frame.render_widget(canvas, canvas_area);

    let mid = (x_bounds.0 + x_bounds.1) / 2.0;
    let axis_label = format!("{:.1}    {:.1}    {:.1} m/z", x_bounds.0, mid, x_bounds.1);
    let axis = Paragraph::new(axis_label)
        .style(Style::default().fg(Color::DarkGray))
        .alignment(Alignment::Center);
    frame.render_widget(axis, axis_area);
}

fn plot_points_cached<'a>(
    app: &'a mut SpectraApp,
    idx: u32,
    data: &SpectrumData<'_>,
    bounds: (f64, f64),
    mode: PlotRenderMode,
    canvas_area: Rect,
) -> (&'a [(f64, f64)], f64) {
    let canvas_width = canvas_area.width;

    let zoom = app.zoom;
    let normalize = app.normalize;
    let hit = matches!(
        app.plot_cache.as_ref(),
        Some(cache)
            if cache.idx == idx
                && cache.zoom == zoom
                && cache.normalize == normalize
                && cache.mode == mode
                && cache.canvas_width == canvas_width
                && approx_bounds_eq(cache.x_bounds, bounds)
    );
    if hit {
        let cache = app
            .plot_cache
            .as_ref()
            .expect("plot_cache checked as Some above");
        return (&cache.points, cache.y_max);
    }

    let bins = plot_bins_for_canvas(canvas_area, mode);
    let (points, y_max) = downsample_max_per_bin(data, bounds, app.normalize, bins);

    if spectra_debug_target(idx) {
        spectra_debug_log(format!(
            "plot idx={idx} mode={} zoom={} normalize={} bounds=[{:.4},{:.4}] bins={bins} points={} y_max={:.4}",
            mode.label(),
            zoom.label(),
            normalize,
            bounds.0,
            bounds.1,
            points.len(),
            y_max
        ));
        for (i, (x, y)) in points.iter().take(10).enumerate() {
            spectra_debug_log(format!("plot sample[{i}] x={x:.4} y={y:.6}"));
        }
    }

    app.plot_cache = Some(PlotCache {
        idx,
        zoom: app.zoom,
        normalize: app.normalize,
        mode,
        canvas_width,
        x_bounds: bounds,
        points,
        y_max,
    });

    let cache = app
        .plot_cache
        .as_ref()
        .expect("plot_cache set immediately above");
    (&cache.points, cache.y_max)
}

fn downsample_max_per_bin(
    data: &SpectrumData<'_>,
    bounds: (f64, f64),
    normalize: bool,
    bins: usize,
) -> (Vec<(f64, f64)>, f64) {
    let mut min_x = bounds.0;
    let mut max_x = bounds.1;
    if min_x > max_x {
        std::mem::swap(&mut min_x, &mut max_x);
    }

    let span = (max_x - min_x).max(1e-9);
    let bins = bins.clamp(16, MAX_PLOT_BINS);

    let mut best_y: Vec<f64> = vec![f64::NEG_INFINITY; bins];
    let mut best_x: Vec<f64> = vec![0.0; bins];
    let mut has: Vec<bool> = vec![false; bins];

    for (&mz, &inten) in data.mz.iter().zip(data.intensity.iter()) {
        if mz < min_x || mz > max_x {
            continue;
        }
        let inten = inten as f64;
        if !inten.is_finite() {
            continue;
        }

        let frac = ((mz - min_x) / span).clamp(0.0, 1.0);
        let mut bin = (frac * (bins as f64)) as usize;
        if bin >= bins {
            bin = bins - 1;
        }

        if !has[bin] || inten > best_y[bin] {
            has[bin] = true;
            best_y[bin] = inten;
            best_x[bin] = mz;
        }
    }

    let mut pts: Vec<(f64, f64)> = has
        .iter()
        .enumerate()
        .filter_map(|(i, present)| {
            if *present {
                Some((best_x[i], best_y[i]))
            } else {
                None
            }
        })
        .collect();

    if pts.is_empty() {
        return (pts, 1.0);
    }

    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    let max_intensity: f64 = pts
        .iter()
        .fold(0.0_f64, |acc, (_, y)| acc.max(*y))
        .max(1e-6_f64);

    if normalize && max_intensity > 0.0 {
        for (_, y) in pts.iter_mut() {
            *y /= max_intensity;
        }
        (pts, 1.1)
    } else {
        (pts, max_intensity * 1.1)
    }
}

fn chromatogram_cached<'a>(app: &'a mut SpectraApp, canvas_width: u16) -> Option<&'a ChromCache> {
    let normalize = app.normalize;
    let hit = matches!(
        app.chrom_cache.as_ref(),
        Some(cache) if cache.normalize == normalize && cache.canvas_width == canvas_width
    );
    if hit {
        return app.chrom_cache.as_ref();
    }

    let mut rt_min = f64::INFINITY;
    let mut rt_max = -f64::INFINITY;
    for meta in app.metas.iter().filter(|m| m.ms_level <= 1) {
        let Some(rt) = meta.rt_minutes else {
            continue;
        };
        let rt = rt as f64;
        rt_min = rt_min.min(rt);
        rt_max = rt_max.max(rt);
    }

    if !rt_min.is_finite() || !rt_max.is_finite() {
        app.chrom_cache = None;
        return None;
    }

    if rt_min >= rt_max {
        rt_max = rt_min + 1.0;
    }
    let rt_bounds = (rt_min, rt_max);

    let bins = (canvas_width as usize)
        .saturating_mul(4)
        .max(256)
        .clamp(16, MAX_PLOT_BINS);

    let tic_samples = app
        .metas
        .iter()
        .filter(|m| m.ms_level <= 1)
        .filter_map(|m| Some((m.rt_minutes? as f64, m.tic? as f64)));
    let (tic_points, tic_y_max) =
        downsample_points_max_per_bin(tic_samples, rt_bounds, normalize, bins);

    let bpc_samples = app
        .metas
        .iter()
        .filter(|m| m.ms_level <= 1)
        .filter_map(|m| Some((m.rt_minutes? as f64, m.base_peak_intensity? as f64)));
    let (bpc_points, bpc_y_max) =
        downsample_points_max_per_bin(bpc_samples, rt_bounds, normalize, bins);

    app.chrom_cache = Some(ChromCache {
        normalize,
        canvas_width,
        rt_bounds,
        tic_points,
        tic_y_max,
        bpc_points,
        bpc_y_max,
    });

    app.chrom_cache.as_ref()
}

fn base_peak_map_cached<'a>(app: &'a mut SpectraApp, canvas_area: Rect) -> Option<&'a MapCache> {
    let normalize = app.normalize;
    let canvas_size = (canvas_area.width, canvas_area.height);

    let hit = matches!(
        app.map_cache.as_ref(),
        Some(cache) if cache.normalize == normalize && cache.canvas_size == canvas_size
    );
    if hit {
        return app.map_cache.as_ref();
    }

    let mut rt_min = f64::INFINITY;
    let mut rt_max = -f64::INFINITY;
    let mut mz_min = f64::INFINITY;
    let mut mz_max = -f64::INFINITY;

    let mut max_intensity = 0.0_f64;
    let mut points = 0usize;

    for meta in app.metas.iter().filter(|m| m.ms_level <= 1) {
        let (Some(rt), Some(mz), Some(intensity)) =
            (meta.rt_minutes, meta.base_peak_mz, meta.base_peak_intensity)
        else {
            continue;
        };
        let rt = rt as f64;
        let mz = mz as f64;
        let intensity = intensity as f64;
        if !rt.is_finite() || !mz.is_finite() || !intensity.is_finite() || intensity <= 0.0 {
            continue;
        }

        points = points.saturating_add(1);
        rt_min = rt_min.min(rt);
        rt_max = rt_max.max(rt);
        mz_min = mz_min.min(mz);
        mz_max = mz_max.max(mz);
        max_intensity = max_intensity.max(intensity);
    }

    if points == 0
        || !rt_min.is_finite()
        || !rt_max.is_finite()
        || !mz_min.is_finite()
        || !mz_max.is_finite()
    {
        app.map_cache = None;
        return None;
    }

    if rt_min >= rt_max {
        rt_max = rt_min + 1.0;
    }
    if mz_min >= mz_max {
        mz_max = mz_min + 1.0;
    }

    let rt_bounds = (rt_min, rt_max);
    let mz_bounds = (mz_min, mz_max);

    let x_bins = (canvas_area.width.max(1) as usize)
        .saturating_mul(2)
        .clamp(16, 1500);
    let y_bins = (canvas_area.height.max(1) as usize)
        .saturating_mul(4)
        .clamp(16, 1500);
    let cell_count = x_bins.saturating_mul(y_bins).max(1);

    let mut best_i = vec![f64::NEG_INFINITY; cell_count];
    let mut best_rt = vec![0.0_f64; cell_count];
    let mut best_mz = vec![0.0_f64; cell_count];
    let mut has = vec![false; cell_count];

    let rt_span = (rt_bounds.1 - rt_bounds.0).max(1e-9);
    let mz_span = (mz_bounds.1 - mz_bounds.0).max(1e-9);

    for meta in app.metas.iter().filter(|m| m.ms_level <= 1) {
        let (Some(rt), Some(mz), Some(intensity)) =
            (meta.rt_minutes, meta.base_peak_mz, meta.base_peak_intensity)
        else {
            continue;
        };
        let rt = rt as f64;
        let mz = mz as f64;
        let intensity = intensity as f64;
        if !rt.is_finite() || !mz.is_finite() || !intensity.is_finite() || intensity <= 0.0 {
            continue;
        }

        let x_frac = ((rt - rt_bounds.0) / rt_span).clamp(0.0, 1.0);
        let y_frac = ((mz - mz_bounds.0) / mz_span).clamp(0.0, 1.0);

        let mut xi = (x_frac * (x_bins as f64)) as usize;
        if xi >= x_bins {
            xi = x_bins - 1;
        }
        let mut yi = (y_frac * (y_bins as f64)) as usize;
        if yi >= y_bins {
            yi = y_bins - 1;
        }

        let idx = yi.saturating_mul(x_bins).saturating_add(xi);
        if idx >= cell_count {
            continue;
        }
        if !has[idx] || intensity > best_i[idx] {
            has[idx] = true;
            best_i[idx] = intensity;
            best_rt[idx] = rt;
            best_mz[idx] = mz;
        }
    }

    let mut buckets: [Vec<(f64, f64)>; 5] = std::array::from_fn(|_| Vec::new());
    let denom = if normalize {
        (1.0 + max_intensity).ln().max(1e-12)
    } else {
        max_intensity.max(1e-12)
    };

    for i in 0..cell_count {
        if !has[i] {
            continue;
        }
        let intensity = best_i[i];
        if !intensity.is_finite() || intensity <= 0.0 {
            continue;
        }

        let frac = if normalize {
            (1.0 + intensity).ln() / denom
        } else {
            (intensity / denom).clamp(0.0, 1.0)
        };
        let mut level = (frac * 5.0).floor() as usize;
        if level >= 5 {
            level = 4;
        }
        buckets[level].push((best_rt[i], best_mz[i]));
    }

    app.map_cache = Some(MapCache {
        normalize,
        canvas_size,
        rt_bounds,
        mz_bounds,
        max_intensity,
        buckets,
    });

    app.map_cache.as_ref()
}

fn downsample_points_max_per_bin(
    samples: impl IntoIterator<Item = (f64, f64)>,
    bounds: (f64, f64),
    normalize: bool,
    bins: usize,
) -> (Vec<(f64, f64)>, f64) {
    let mut min_x = bounds.0;
    let mut max_x = bounds.1;
    if min_x > max_x {
        std::mem::swap(&mut min_x, &mut max_x);
    }

    let span = (max_x - min_x).max(1e-9);
    let bins = bins.clamp(16, MAX_PLOT_BINS);

    let mut best_y: Vec<f64> = vec![f64::NEG_INFINITY; bins];
    let mut best_x: Vec<f64> = vec![0.0; bins];
    let mut has: Vec<bool> = vec![false; bins];

    for (x, y) in samples {
        if x < min_x || x > max_x || !x.is_finite() || !y.is_finite() {
            continue;
        }

        let frac = ((x - min_x) / span).clamp(0.0, 1.0);
        let mut bin = (frac * (bins as f64)) as usize;
        if bin >= bins {
            bin = bins - 1;
        }

        if !has[bin] || y > best_y[bin] {
            has[bin] = true;
            best_y[bin] = y;
            best_x[bin] = x;
        }
    }

    let mut pts: Vec<(f64, f64)> = has
        .iter()
        .enumerate()
        .filter_map(|(i, present)| {
            if *present {
                Some((best_x[i], best_y[i]))
            } else {
                None
            }
        })
        .collect();

    if pts.is_empty() {
        return (pts, 1.0);
    }

    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    let max_y: f64 = pts
        .iter()
        .fold(0.0_f64, |acc, (_, y)| acc.max(*y))
        .max(1e-6_f64);

    if normalize && max_y > 0.0 {
        for (_, y) in pts.iter_mut() {
            *y /= max_y;
        }
        (pts, 1.1)
    } else {
        (pts, max_y * 1.1)
    }
}

fn compute_bounds(zoom: Zoom, meta: &SpectrumMeta, data: &SpectrumData<'_>) -> (f64, f64) {
    let min_mz = meta.mz_min.unwrap_or_else(|| {
        data.mz
            .iter()
            .copied()
            .fold(f64::INFINITY, |acc, v| acc.min(v))
    });
    let max_mz = meta.mz_max.unwrap_or_else(|| {
        data.mz
            .iter()
            .copied()
            .fold(-f64::INFINITY, |acc, v| acc.max(v))
    });
    let full = if min_mz.is_finite() && max_mz.is_finite() && min_mz < max_mz {
        (min_mz, max_mz)
    } else {
        (0.0, 2000.0)
    };

    match zoom {
        Zoom::Full => full,
        Zoom::Low => (0.0, 500.0),
        Zoom::Mid => (500.0, 1000.0),
        Zoom::High => (1000.0, 2000.0),
        Zoom::Precursor => {
            if let Some(pmz) = meta.precursor_mz {
                let center = pmz;
                (center - 75.0, center + 75.0)
            } else {
                full
            }
        }
    }
}

fn resolve_plot_mode(
    mode: PlotMode,
    continuity: SignalContinuity,
    points: usize,
) -> PlotRenderMode {
    match mode {
        PlotMode::Auto => match continuity {
            SignalContinuity::Profile => PlotRenderMode::Line,
            SignalContinuity::Unknown if points > 5_000 => PlotRenderMode::Line,
            _ => PlotRenderMode::Sticks,
        },
        PlotMode::Sticks => PlotRenderMode::Sticks,
        PlotMode::Line => PlotRenderMode::Line,
    }
}

fn plot_bins_for_canvas(canvas_area: Rect, mode: PlotRenderMode) -> usize {
    let width = canvas_area.width.max(1) as usize;
    let base = match mode {
        PlotRenderMode::Sticks => width.saturating_mul(4).max(128),
        PlotRenderMode::Line => width.saturating_mul(8).max(256),
    };
    base.clamp(16, MAX_DRAW_PEAKS)
}

fn approx_bounds_eq(a: (f64, f64), b: (f64, f64)) -> bool {
    approx_eq(a.0, b.0) && approx_eq(a.1, b.1)
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9
}

fn continuity_short(continuity: SignalContinuity) -> &'static str {
    match continuity {
        SignalContinuity::Centroid => "C",
        SignalContinuity::Profile => "P",
        SignalContinuity::Unknown => "?",
    }
}

fn export_current_spectrum_svg(app: &mut SpectraApp) -> anyhow::Result<PathBuf> {
    let meta = app
        .current_meta()
        .ok_or_else(|| anyhow::anyhow!("no spectrum selected"))?;

    let data = app
        .data(meta.idx)
        .ok_or_else(|| anyhow::anyhow!("spectrum not loaded (press Enter first)"))?;

    let x_bounds = compute_bounds(app.zoom, meta, &data);
    let mode = resolve_plot_mode(app.plot_mode, meta.continuity, data.mz.len());
    let (points, y_max) = downsample_max_per_bin(&data, x_bounds, app.normalize, SVG_BINS);

    fs::create_dir_all("exports").context("failed to create exports/")?;

    let scan_id = sanitize_filename_component(&meta.scan_id);
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    let filename = format!(
        "spectrum_{:05}_{}_ms{}_{}_{}_{}.svg",
        meta.idx,
        scan_id,
        meta.ms_level,
        app.zoom.label(),
        mode.label(),
        ts
    );
    let out_path = PathBuf::from("exports").join(filename);

    write_spectrum_svg(
        &out_path,
        SVG_WIDTH,
        SVG_HEIGHT,
        meta,
        mode,
        x_bounds,
        y_max,
        &points,
        app.normalize,
    )?;

    Ok(out_path)
}

fn export_current_spectrum_png(app: &mut SpectraApp) -> anyhow::Result<PathBuf> {
    use std::io::ErrorKind;
    use std::process::{Command, Stdio};

    let svg_path = export_current_spectrum_svg(app)?;
    let png_path = svg_path.with_extension("png");

    let inkscape = Command::new("inkscape")
        .arg(&svg_path)
        .arg("--export-type=png")
        .arg("--export-filename")
        .arg(&png_path)
        .arg("--export-area-drawing")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
    match inkscape {
        Ok(status) if status.success() => return Ok(png_path),
        Ok(_) => {}
        Err(err) if err.kind() == ErrorKind::NotFound => {}
        Err(err) => return Err(err).context("failed to run inkscape for PNG export"),
    }

    let convert = Command::new("convert")
        .arg(&svg_path)
        .arg(&png_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
    match convert {
        Ok(status) if status.success() => Ok(png_path),
        Ok(_) => Err(anyhow::anyhow!(
            "PNG conversion failed via inkscape/convert (SVG saved at {})",
            svg_path.display()
        )),
        Err(err) if err.kind() == ErrorKind::NotFound => Err(anyhow::anyhow!(
            "no SVG→PNG converter found (tried inkscape/convert); SVG saved at {}",
            svg_path.display()
        )),
        Err(err) => Err(err).context("failed to run convert for PNG export"),
    }
}

fn sanitize_filename_component(input: &str) -> String {
    let mut out = String::with_capacity(input.len().min(80));
    for ch in input.chars().take(80) {
        match ch {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '-' | '_' | '.' => out.push(ch),
            _ => out.push('_'),
        }
    }
    if out.is_empty() {
        "scan".to_string()
    } else {
        out
    }
}

fn write_spectrum_svg(
    path: &Path,
    width: u32,
    height: u32,
    meta: &SpectrumMeta,
    mode: PlotRenderMode,
    x_bounds: (f64, f64),
    y_max: f64,
    points: &[(f64, f64)],
    normalize: bool,
) -> anyhow::Result<()> {
    let mut file =
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?;

    let margin_left = 70.0;
    let margin_right = 20.0;
    let margin_top = 30.0;
    let margin_bottom = 60.0;
    let w = width as f64;
    let h = height as f64;
    let plot_w = (w - margin_left - margin_right).max(1.0);
    let plot_h = (h - margin_top - margin_bottom).max(1.0);

    let x0 = x_bounds.0.min(x_bounds.1);
    let x1 = x_bounds.0.max(x_bounds.1);
    let x_span = (x1 - x0).max(1e-9);
    let y_span = y_max.max(1e-9);

    let title = format!(
        "Scan {} | ms{} | rt={}min | {} | {}",
        meta.scan_id,
        meta.ms_level,
        meta.rt_minutes
            .map(|v| format!("{v:.2}"))
            .unwrap_or_else(|| "-".to_string()),
        meta.continuity,
        if normalize { "normalized" } else { "raw" }
    );

    writeln!(
        file,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"##
    )?;
    writeln!(
        file,
        r##"<rect x="0" y="0" width="{width}" height="{height}" fill="white"/>"##
    )?;
    writeln!(
        file,
        r##"<text x="{x}" y="{y}" font-family="Helvetica, Arial, sans-serif" font-size="16" fill="#111">{text}</text>"##,
        x = margin_left,
        y = 20,
        text = escape_xml(&title),
    )?;

    // Plot frame.
    writeln!(
        file,
        r##"<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="none" stroke="#333" stroke-width="1"/>"##,
        x = margin_left,
        y = margin_top,
        w = plot_w,
        h = plot_h,
    )?;

    let to_px = |x: f64, y: f64| -> (f64, f64) {
        let xf = ((x - x0) / x_span).clamp(0.0, 1.0);
        let yf = (y / y_span).clamp(0.0, 1.0);
        let px = margin_left + xf * plot_w;
        let py = margin_top + (1.0 - yf) * plot_h;
        (px, py)
    };

    // Series.
    match mode {
        PlotRenderMode::Sticks => {
            for &(x, y) in points.iter() {
                let (px, py) = to_px(x, y);
                let (px0, py0) = to_px(x, 0.0);
                writeln!(
                    file,
                    r##"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="#00a8c6" stroke-width="1" />"##,
                    x1 = px0,
                    y1 = py0,
                    x2 = px,
                    y2 = py,
                )?;
            }
        }
        PlotRenderMode::Line => {
            let mut d = String::new();
            for (i, &(x, y)) in points.iter().enumerate() {
                let (px, py) = to_px(x, y);
                if i == 0 {
                    d.push_str(&format!("M{px:.2},{py:.2}"));
                } else {
                    d.push_str(&format!(" L{px:.2},{py:.2}"));
                }
            }
            writeln!(
                file,
                r##"<path d="{d}" fill="none" stroke="#00a8c6" stroke-width="1"/>"##,
                d = d
            )?;
        }
    }

    // Axes labels.
    let (x_left, _) = to_px(x0, 0.0);
    let (x_mid, _) = to_px((x0 + x1) / 2.0, 0.0);
    let (x_right, _) = to_px(x1, 0.0);
    let y_label = margin_top + plot_h + 35.0;

    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#333" text-anchor="middle">{text}</text>"##,
        x = x_left,
        y = y_label,
        text = format!("{x0:.2}"),
    )?;
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#333" text-anchor="middle">{text}</text>"##,
        x = x_mid,
        y = y_label,
        text = format!("{:.2}", (x0 + x1) / 2.0),
    )?;
    writeln!(
        file,
        r##"<text x="{x:.2}" y="{y:.2}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#333" text-anchor="middle">{text}</text>"##,
        x = x_right,
        y = y_label,
        text = format!("{x1:.2} m/z"),
    )?;

    writeln!(file, "</svg>")?;
    Ok(())
}

fn escape_xml(input: &str) -> String {
    input
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('\"', "&quot;")
        .replace('\'', "&apos;")
}

fn init_terminal() -> anyhow::Result<CrosstermTerminal> {
    enable_raw_mode().context("failed to enable raw mode")?;
    let mut stdout = std::io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)
        .context("failed to enter alternate screen")?;
    let backend = CrosstermBackend::new(stdout);
    Terminal::new(backend).context("failed to build terminal")
}

fn restore_terminal(terminal: &mut CrosstermTerminal) -> anyhow::Result<()> {
    disable_raw_mode().context("failed to disable raw mode")?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )
    .context("failed to leave alternate screen")?;
    terminal.show_cursor().context("failed to show cursor")
}
