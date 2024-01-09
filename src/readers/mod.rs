#![doc(hidden)]
//! Library of readers and common functions for various file formats

// internal modules
use crate::mesh::Mesh;
use crate::point::Point;
use crate::posvol::Posvol;

// standard library
use std::path::Path;

// external crates
use anyhow::Result;

// files under the readers module
mod meshtal_file;
pub mod parsers;
mod points_file;
mod posvol_file;

// inline important the mesh-related modules for a nice API
#[doc(inline)]
pub use crate::readers::meshtal_file::MeshtalReader;

#[doc(inline)]
pub use crate::readers::points_file::PointsFileReader;

#[doc(inline)]
pub use crate::readers::posvol_file::PosvolReader;

/// Read all meshes in a meshtal file
///
/// Returns a result containing a vector of [Mesh] structs extracted from the
/// file at `path` by the parser.
///
/// - `path` - Path to the meshtal file, can be [&str], [String], [Path], etc...
///
/// Example
/// ```ignore
/// // Read every mesh contained in the file
/// let mesh_tallies: Vec<Mesh> = meshtal::read_meshtal("path/to/meshtal.msht")?;
/// ```
pub fn read_meshtal<P: AsRef<Path>>(path: P) -> Result<Vec<Mesh>> {
    let path: &Path = Path::new(path.as_ref());
    let mut reader = MeshtalReader::new();
    reader.disable_progress();
    reader.parse(path)
}

/// Read only the specified mesh from a meshtal file
///
/// Returns a result of the targeted [Mesh] if it was successfully
/// extracted from the file at `path`.
///
/// - `path` - Path to the meshtal file, can be [&str], [String], [Path], etc...
/// - `target` - Tally number of interest
///
/// Example
/// ```ignore
/// // Read only tally 104 (i.e. FMESH104) from the file
/// let mesh: Mesh = meshtal::read_meshtal_target("path/to/meshtal.msht", 104)?;
/// ```
pub fn read_meshtal_target<P: AsRef<Path>>(path: P, target: u32) -> Result<Mesh> {
    let path: &Path = Path::new(path.as_ref());
    let mut reader = MeshtalReader::new();
    reader.disable_progress();
    reader.set_target_id(target);
    let mut mesh_list = reader.parse(path)?;
    Ok(mesh_list.remove(0))
}

/// Deserialise binary posvol file
///
/// Returns a Result containing a Posvol struct with all the information
/// extracted from a CuV posvol file at `path``.
pub fn read_posvol_file<P: AsRef<Path>>(path: P) -> Result<Posvol> {
    let path: &Path = Path::new(path.as_ref());
    let reader = PosvolReader::new();
    reader.parse(path)
}

/// Read points file for the pointextract tool
///
/// Returns a result containing a vector of [Point]
/// structs extracted from the provided input file.
///
/// - `path` - Path to a points file, can be [&str], [String], [Path], etc...
///
/// The input file (default `points.txt`) is interpreted with the following
/// rules for a line:
///
/// | Example line               | Interpretation           |
/// | -------------------------- | ------------------------ |
/// | Starts with `#`            | comment                  |
/// | `rzt`, `cyl`, `xyz`, `rec` | geometry keyword         |
/// | 1.0 2.0 3.0                | i, j, k                  |
/// | 1e2  1.0 2.0 3.0           | energy, i, j, k          |
/// | 1e2 total 1.0 2.0 3.0      | energy, time, i, j, k    |
///
/// Anything else is ignored. For examples see `data/points.txt`.
///
/// Example
/// ```ignore
/// // Read any valid user points from a text file
/// let points: Vec<Point> = meshtal::read_points_file("path/to/points.txt")?;
/// ```
pub fn read_points_file<P: AsRef<Path>>(path: P) -> Result<Vec<Point>> {
    let path: &Path = Path::new(path.as_ref());
    let mut reader = PointsFileReader::new();
    reader.parse(path)
}
