#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct CoordinateRange<T> {
    pub(crate) start: T,
    pub(crate) end: T,
}

impl CoordinateRange<f64> {
    pub(crate) fn new(start: f64, end: f64) -> Self {
        Self { start, end }
    }

    pub(crate) fn min(self) -> f64 {
        self.start.min(self.end)
    }

    pub(crate) fn max(self) -> f64 {
        self.start.max(self.end)
    }

    pub(crate) fn span(self) -> f64 {
        self.end - self.start
    }

    pub(crate) fn normalized_position(self, value: f64) -> f64 {
        let span = self.span();
        if !span.is_finite() || span.abs() <= f64::EPSILON {
            return 0.0;
        }
        (value - self.start) / span
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Scale {
    pub(crate) domain: CoordinateRange<f64>,
    pub(crate) range: CoordinateRange<f64>,
}

impl Scale {
    pub(crate) fn new(domain: CoordinateRange<f64>, range: CoordinateRange<f64>) -> Self {
        Self { domain, range }
    }

    pub(crate) fn transform(self, value: f64) -> f64 {
        let position = self.domain.normalized_position(value);
        self.range.start + position * self.range.span()
    }
}

#[cfg(test)]
mod tests {
    use super::{CoordinateRange, Scale};

    #[test]
    fn coordinate_range_handles_descending_ranges() {
        let range = CoordinateRange::new(100.0, 0.0);
        assert!((range.normalized_position(20.0) - 0.8).abs() < 1e-9);
    }

    #[test]
    fn scale_transforms_into_target_range() {
        let scale = Scale::new(
            CoordinateRange::new(10.0, 20.0),
            CoordinateRange::new(100.0, 200.0),
        );
        assert!((scale.transform(15.0) - 150.0).abs() < 1e-9);
    }
}
