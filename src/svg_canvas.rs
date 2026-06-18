use crate::scale::{CoordinateRange, Scale};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AxisOrientation {
    Bottom,
    Left,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AxisTickLabelStyle {
    Precision(usize),
    Scientific(usize),
}

impl AxisTickLabelStyle {
    pub(crate) fn format(self, value: f64) -> String {
        match self {
            Self::Precision(precision) => format!("{value:.precision$}"),
            Self::Scientific(precision) => format!("{value:.precision$e}"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct AxisProps {
    orientation: AxisOrientation,
    label: &'static str,
    tick_label_style: AxisTickLabelStyle,
}

impl AxisProps {
    pub(crate) fn new(orientation: AxisOrientation, label: &'static str) -> Self {
        Self {
            orientation,
            label,
            tick_label_style: AxisTickLabelStyle::Precision(2),
        }
    }

    pub(crate) fn with_tick_label_style(mut self, tick_label_style: AxisTickLabelStyle) -> Self {
        self.tick_label_style = tick_label_style;
        self
    }

    pub(crate) fn format_tick(self, value: f64) -> String {
        self.tick_label_style.format(value)
    }

    pub(crate) fn label(self) -> &'static str {
        self.label
    }

    #[allow(dead_code)]
    pub(crate) fn orientation(self) -> AxisOrientation {
        self.orientation
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SvgCanvas {
    left: f64,
    top: f64,
    width: f64,
    height: f64,
    x_scale: Scale,
    y_scale: Scale,
}

impl SvgCanvas {
    pub(crate) fn new(
        left: f64,
        top: f64,
        width: f64,
        height: f64,
        x_domain: CoordinateRange<f64>,
        y_domain: CoordinateRange<f64>,
    ) -> Self {
        let x_scale = Scale::new(x_domain, CoordinateRange::new(left, left + width));
        let y_scale = Scale::new(y_domain, CoordinateRange::new(top + height, top));
        Self {
            left,
            top,
            width,
            height,
            x_scale,
            y_scale,
        }
    }

    pub(crate) fn left(self) -> f64 {
        self.left
    }

    pub(crate) fn top(self) -> f64 {
        self.top
    }

    pub(crate) fn width(self) -> f64 {
        self.width
    }

    pub(crate) fn height(self) -> f64 {
        self.height
    }

    pub(crate) fn right(self) -> f64 {
        self.left + self.width
    }

    pub(crate) fn bottom(self) -> f64 {
        self.top + self.height
    }

    pub(crate) fn x(self, value: f64) -> f64 {
        self.x_scale.transform(value)
    }

    pub(crate) fn y(self, value: f64) -> f64 {
        self.y_scale.transform(value)
    }

    pub(crate) fn transform(self, x: f64, y: f64) -> (f64, f64) {
        (self.x(x), self.y(y))
    }
}

#[cfg(test)]
mod tests {
    use super::{AxisOrientation, AxisProps, AxisTickLabelStyle, SvgCanvas};
    use crate::scale::CoordinateRange;

    #[test]
    fn svg_canvas_transforms_points_into_plot_space() {
        let canvas = SvgCanvas::new(
            10.0,
            20.0,
            200.0,
            100.0,
            CoordinateRange::new(0.0, 100.0),
            CoordinateRange::new(0.0, 10.0),
        );
        let (x, y) = canvas.transform(25.0, 5.0);
        assert!((x - 60.0).abs() < 1e-9);
        assert!((y - 70.0).abs() < 1e-9);
    }

    #[test]
    fn axis_props_format_ticks() {
        let axis = AxisProps::new(AxisOrientation::Left, "Intensity")
            .with_tick_label_style(AxisTickLabelStyle::Scientific(2));
        assert_eq!(axis.format_tick(1234.0), "1.23e3");
    }
}
