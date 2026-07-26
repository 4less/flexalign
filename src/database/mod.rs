pub mod common;
// `build` held the legacy hash-index builder (flexmap::build::hash / FlexmapHash); the live path
// builds via `flexmap::build::default_build` in `flexmap::DB`. Disabled after flexmap dropped the
// hash API. The module file is retained but not compiled.
// pub mod build;
pub mod flexmap;
