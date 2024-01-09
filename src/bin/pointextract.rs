//! Command line extraction of point data from meshes
//!
//! Provides a quick method for finding data at specific points in a 3D mesh
//! with support for all mesh types and geometries. For consistency with legacy
//! tools and large numbers of points, an input file option is also supported.
//!
//! # Usage
//!
//! ```text
//! Usage: pointextract <meshtal> <number> [options]
//! ```
//!
//! Help is printed with the `-h` flag, and `--help` will show examples, default
//! values, examples, and any important behaviour.
//!
//! ### Definition of a 'point'
//!
//! A `Point` in may be defined three ways. Using cartesian `xyz` as an example:
//!
//! - (x, y, z)
//! - (energy, x, y, z)
//! - (energy, time, x, y, z)
//!
//! where are `energy`/`time` groups are either a value or 'total', and `x`,
//! `y`, `z` are the location to search for.
//!
//! The 'Total' group is used by default for omitted `energy` and `time` groups.
//!
//! Note that no matter what coordinates are provided, they are converted into
//! the appropriate coordinate system in the background. For example, `xyz`
//! coordinates are converted to `rzt` if the target mesh is cylindrical.
//!
//! Points can be retrieved in two ways:
//!
//! - Single point via the `-p`/`--point` argument
//! - Multiple points via the `-f`/`--file` argument
//!
//! ### Single point (--point / --type)
//!
//! To quickly get a single point:
//!
//! ```bash
//! # Find the point at (x,y,z) = (1.0, 2.0, 3.0)
//! pointextract /path/to/meshtal.msht 104 -p 1.0 2.0 3.0
//!
//! # Find the same point, but for the 100.0 MeV energy group
//! pointextract /path/to/meshtal.msht 104 -p 1.0E+02 1.0 2.0 3.0
//!
//! # Find the same point, but for the 'Total' energy group and 6.0e+15 time group
//! pointextract /path/to/meshtal.msht 104 -p total 6e+15 1.0 2.0 3.0
//! ```
//!
//! This works for both the `Rectangular` and `Cylindrical` mesh types.
//!
//! Coordinates default to being XYZ cartesian, but the `--type` argument allows
//! this to be specified explicitly
//!
//! ```bash
//! # Find the point at (r,z,t) = (1.0, 2.0, 90.0)
//! pointextract /path/to/meshtal.msht 104 -p 1.0 2.0 0.53 --type rzt
//! ```
//!
//! ### Multiple points (--file)
//!
//! The input file (default `points.txt`) is interpreted with the following
//! rules for a line:
//!
//! | Example line               | Interpretation           |
//! | -------------------------- | ------------------------ |
//! | Starts with `#`            | comment                  |
//! | `rzt`, `cyl`, `xyz`, `rec` | geometry keyword         |
//! | 1.0 2.0 3.0                | i, j, k                  |
//! | 1e2  1.0 2.0 3.0           | energy, i, j, k          |
//! | 1e2 total 1.0 2.0 3.0      | energy, time, i, j, k    |
//!
//! Anything else is ignored. For an example file see `data/points.txt`, though
//! a simplified input is shown below.
//!
//! ```bash
//! # this is a line comment in points.txt
//!
//! xyz                         # points below explicitly interpreted as cartesian  
//! 1.0 5.0 7.0                 # 'Total' energy, 'Total' time, (x, y, z)
//! total 1.0 5.0 7.0           # 'Total' energy, 'Total' time, (x, y, z)
//! total total 1.0 5.0 7.0     # 'Total' energy, 'Total' time, (x, y, z)
//!
//! rzt                         # points below explicitly interpreted as cylindrical  
//! 4.0 1.0 5.0 0.5             #  4 MeV  energy, 'Total' time, (r, z, t)
//! 4.0 1e16 1.0 5.0 0.5        #  4 MeV  energy,  1e16   time, (r, z, t)
//! ```
//!
//! This is used as
//!
//! ```bash
//! # Find all points in file
//! pointextract /path/to/meshtal.msht 104 --file points.txt
//! ```
//!
//! It is fine to mix and match coordinates in the same file because all points
//! are converted into the appropriate coordinate system in the background.
//!
//! Lines that can not be parsed into a `Point` are ignored, and warnings are
//! raised for invalid points that are outside of the mesh bounds.
//!
//! ### Result outputs
//!
//! Results are written to `results.dat` by default but this can be renamed as
//! needed.
//!
//! ```bash
//! # Search for all locations in points.txt, output results to 'myoutput.txt'
//! pointextract /path/to/meshtal.msht 104 --file points.txt --output myoutput.txt
//! ```
//!
//! These may also be written to the terminal directly with the `-d`/`--dump`
//! flag.
//!
//! Both the parsed user input and search results are presented in tables. For
//! example, a 'points.txt' file may have the following
//!
//! ```bash
//! total 1.111e16 0.5 0.0 -1.9    # (energy, time, x, y, z), inside of mesh bounds
//! total 1.111e16 50.5 0.0 -1.9   # (energy, time, x, y, z), outside of mesh bounds
//! ```
//!
//! Tables show the user provided points to search for, and detail of the voxel containing these points.
//!
//! ```text
//!                             Points to search
//! id     energy        time        i_coord      j_coord      k_coord   system
//! -----------------------------------------------------------------------------
//! 0      Total     1.11100e+16  5.00000e-01  0.00000e+00 -1.90000e+00   xyz
//! 1      Total     1.11100e+16  5.05000e+01  0.00000e+00 -1.90000e+00   xyz
//!
//!                     Voxels found (fmesh334, Rectangular)
//! id     energy        time     i_coord   j_coord   k_coord    result    error
//! --------------------------------------------------------------------------------
//! 0      Total     1.11100e+16   0.500     0.000    -2.625  1.37928e-02 0.0092
//! 1                             Not found in mesh
//! ```
//!
//! Note the second point was outside of the mesh and failed, which will warn
//! you with a `Not found in mesh` entry. Corresponding rows are numbered for
//! convenience.
//! *Note - long rows are really inconvenient for reading results on a lot of
//! screens, so the choice was made to split up the two tables. A case could be
//! made for doing somthing else with the results for easier parsing.*
//!

// standard library
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

// Crate modules
use meshtal::mesh::{Geometry, Mesh, Voxel};
use meshtal::point::{self, Point};
use meshtal::readers::{parsers, MeshtalReader};
use meshtal::utils::*;

// External crates
use anyhow::{anyhow, Ok, Result};
use clap::{arg, Parser};
use log::*;

#[doc(hidden)]
fn main() -> Result<()> {
    // set up the command line interface and match arguments
    let cli: Cli = Cli::parse();

    // set up logging (+2 to make Info the default)
    let verbosity = cli.verbose + 2;
    logging_init(verbosity, cli.quiet);

    // check points from user before bothering with the mesh
    debug!("Collecting requested points");
    let mut points: Vec<Point> = Vec::new();
    match cli.point.is_empty() {
        false => points.push(parse_cli_point(&cli)?),
        true => {
            let path = cli.file.clone().unwrap_or("points.txt".to_string());
            points.extend(point::read_points_file(path)?)
        }
    }
    check_points(&points)?;
    info!("Point conversion successful");

    // Get the mesh tally
    info!("Reading {}", &cli.meshtal);
    let mesh = try_meshtal_read(&cli)?;
    info!("Mesh read successful");

    // find the voxels
    let voxels = point::find_voxels(&mesh, &points);
    if voxels.is_empty() {
        error!("No valid points found in mesh");
        return Err(anyhow!("No valid points found in Fmesh{}", mesh.id));
    }

    // output to a file and log to console if requested
    let description = f!("fmesh:{}, {:?}", mesh.id, mesh.geometry);
    match cli.output {
        Some(path) => results_to_file(&path, &mesh, &points, &voxels, &description)?,
        None => results_to_file("./results.dat", &mesh, &points, &voxels, &description)?,
    }

    if cli.dump {
        results_to_console(&mesh, &points, &voxels, &description)
    }

    Ok(())
}

/// Extract point data from any mesh in a meshtal file
///
/// All meshtal formats are supported for rectangular and cylindrical meshes.
///
/// If no particular energy/time groups are specified the 'Total' will be used.
///
/// If void record is 'off' for UKAEA CuV tallies, void voxels will have a flux
/// of zero.
///
/// Important: Points exactly on a voxel boundary are assumed to be in the
/// lowest bounds. Averaging of voxels in this case is a work in progress.
///
/// Examples
/// --------
///
///  Typical use:
///     $ pointextract run0.msht 104 -f points.txt -o results.dat
///
///  Point from (x,y,z) coordiantes:
///     $ pointextract run0.msht 104 -p 1.0 2.0 3.0
///
///  Point from (r,z,t) coordiantes:
///     $ pointextract run0.msht 104 -p 1.0 2.0 3.0 -t rzt
///
///  Point with specific energy/time groups:
///     $ pointextract run0.msht 104 -p 1e2 total 1.0 2.0 3.0
///
///  Multiple points from input file:
///     $ pointextract run0.msht 104 -f points.txt
///   
///     points.txt interpretation:
///         '#' or non-keyword   = ignored
///         0.0 0.2 50           = i, j, k
///         1e2 15 2.0 3.1       = energy, i, j, k
///         1e2 total 15 2.0 3.1 = energy, time, i, j, k
///         rzt, cyl, xyz, rec   = keywords
///   
///     All points following a keyword are interpreted as the specified
///     geometry. Otherwise coordiantes assumed to match the mesh type.
///
///     Coordinates are converted into the correct coordinate system
///     in the background during the search.
///
#[allow(rustdoc::invalid_rust_codeblocks)]
#[doc(hidden)]
#[derive(Parser, Debug)]
#[command(
    verbatim_doc_comment,
    arg_required_else_help(true),
    before_help(banner()),
    after_help("Typical use: pointextract run0.msht 104 -f points.txt -o results.dat\n\nNOTE: --help shows more detail and examples"),
    term_width(70),
    hide_possible_values(true),
    override_usage("pointextract <meshtal> <number> [options]")
)]
struct Cli {
    // * Positional
    /// Path to input meshtal file
    #[arg(name = "meshtal")]
    meshtal: String,

    /// Mesh tally identifier
    ///
    /// e.g. 104 for FMESH104:n
    #[arg(name = "number")]
    number: u32,

    // * Optional
    /// Quickly find a single point
    ///
    /// Expressed as (e,t,i,j,k), but only the i, j, k coordinates are required.
    ///     
    /// E.g. -p 1.0 2.0 3.0 as a minimum input.
    ///
    /// Energy and time targets may be ommited for just 'total' bins, or
    /// specified to be explicit.
    ///     
    /// E.g. -p total 1e22 1.0 2.0 3.0
    ///
    /// Assumes coordiantes are appropriate for the mesh type provided by
    /// --type, e.g. '--type rzt' for cylindrical coordinates
    #[arg(help_heading("Point options"))]
    #[arg(short, long)]
    #[arg(allow_negative_numbers(true),
    num_args(1..=5))]
    #[clap(required = false)]
    #[arg(value_name = "coords")]
    point: Vec<String>,

    /// Coordinate system for --point (e.g. xyz)
    ///
    /// Specified the coordinate system to interpret the --point argument by.
    ///
    /// Systems available:
    ///     > xyz = Cartesian coordinate system (default)
    ///     > rzt = Cylindrical coordinate system
    #[arg(help_heading("Point options"))]
    #[arg(short, long, value_enum)]
    #[arg(hide_default_value(true))]
    #[arg(default_value_t = Geometry::Rectangular)]
    #[arg(verbatim_doc_comment)]
    #[arg(id = "type")]
    type_: Geometry,

    /// File containing multiple points
    ///
    /// If no points explicitly given using --point, the tool will
    /// automatically search for a file to read points from instead.
    /// Defaults to "points.txt".
    ///
    ///     points.txt interpretation:
    ///         '#' or non-keyword   = ignored
    ///         0.0 0.2 50           = i, j, k
    ///         1e2 15 2.0 3.1       = energy, i, j, k
    ///         1e2 total 15 2.0 3.1 = energy, time, i, j, k
    ///         rzt, cyl, xyz, rec   = keywords
    ///  
    /// Points following a keyword read as the specified geometry.
    /// Otherwise coordiantes assumed to match the mesh type.
    #[arg(help_heading("Point options"))]
    #[arg(short, long)]
    #[arg(value_name = "path")]
    #[arg(verbatim_doc_comment)]
    file: Option<String>,

    /// Write results to a text file
    ///
    /// If provided, the results of all points are written to this file.
    /// Includes any failed searches for user reference.
    #[arg(help_heading("Output options"))]
    #[arg(short, long)]
    #[clap(required = false)]
    #[arg(allow_negative_numbers(true),
    num_args(0..2))]
    #[arg(value_name = "path")]
    output: Option<String>,

    /// Print results to the terminal
    ///
    /// Tables are broken up into search points and found voxels for
    /// readability. Row indexing is provided for convenience.
    #[arg(help_heading("Output options"))]
    #[arg(short, long)]
    #[clap(required = false)]
    dump: bool,

    // * Flags
    /// Verbose logging (-v, -vv)
    ///
    /// If specified, the default log level of INFO is increased to DEBUG (-v)
    /// or TRACE (-vv). Errors and Warnings are always logged unless in quiet
    /// (-q) mode.
    #[arg(short, long)]
    #[arg(action = clap::ArgAction::Count)]
    verbose: u8,

    /// Supress all log output (overrules --verbose)
    #[arg(short, long)]
    quiet: bool,
}

#[doc(hidden)]
fn try_meshtal_read(cli: &Cli) -> Result<Mesh> {
    let path: &Path = Path::new(&cli.meshtal);

    let mut reader = MeshtalReader::new();
    println!("verbosity {}", cli.verbose);
    if cli.quiet || cli.verbose > 1 {
        reader.disable_progress();
    }
    reader.set_target_id(cli.number);

    let mut mesh = reader.parse(path)?;
    Ok(std::mem::take(&mut mesh[0]))
}

#[doc(hidden)]
fn check_points(points: &[Point]) -> Result<()> {
    if points.is_empty() {
        Err(anyhow!("No valid point input found"))
    } else {
        Ok(())
    }
}

#[doc(hidden)]
fn parse_cli_point(cli: &Cli) -> Result<Point> {
    match &mut parsers::points_file_point(&cli.point.join(" ")) {
        nom::IResult::Ok(data) => {
            let (_, point) = data;
            point.coordinate_type = cli.type_;
            Ok(point.to_owned())
        }
        _ => Err(anyhow!(
            "Failed to parse \"{}\" to a Point",
            &cli.point.join(" ")
        )),
    }
}

/// Write all results to a file
///
/// This will include those that failed for the user's reference.
#[doc(hidden)]
fn results_to_file(
    path: &str,
    mesh: &Mesh,
    points: &[Point],
    voxels: &[Option<Voxel>],
    description: &str,
) -> Result<()> {
    info!("Writing results to {}", path);
    let f = File::create(path)?;
    let mut f = BufWriter::new(f);
    f.write_all(points_table(points).as_bytes())?;
    f.write_all(b"\n\n\n")?;
    f.write_all(voxels_table(mesh, voxels, description).as_bytes())?;
    Ok(())
}

/// Log the results to console
///
/// Will be exactly the same as the table that gets dumped to file
#[doc(hidden)]
fn results_to_console(mesh: &Mesh, points: &[Point], voxels: &[Option<Voxel>], description: &str) {
    println!("\n{}", points_table(points));
    println!("\n\n{}", voxels_table(mesh, voxels, description));
}

/// generates a banner for cli tool consistency
#[doc(hidden)]
fn banner() -> String {
    let mut s = f!("{:-<1$}\n", "", 70);
    s += &f!("{:^70}\n", "Meshtal :: PointExtract");
    s += &f!("{:-<1$}", "", 70);
    s
}

#[doc(hidden)]
fn logging_init(verbosity: u8, quiet: bool) {
    stderrlog::new()
        .modules(vec![
            module_path!(),
            "meshtal::point",
            "meshtal::mesh",
            "meshtal::readers::meshtal_file",
        ])
        .quiet(quiet)
        .verbosity(verbosity as usize)
        .show_level(false)
        .color(stderrlog::ColorChoice::Never)
        .timestamp(stderrlog::Timestamp::Off)
        .init()
        .unwrap();
}

#[doc(hidden)]
fn points_table(points: &[Point]) -> String {
    let mut s = target_heading();
    s += &f!("\n{}\n", target_columns());
    s += &"-".repeat(78);
    points
        .iter()
        .enumerate()
        .for_each(|(i, p)| s += &f!("\n{i:^6}{}", target_row(p)));
    s
}

#[doc(hidden)]
fn voxels_table(mesh: &Mesh, voxels: &[Option<Voxel>], description: &str) -> String {
    let mut s = voxel_heading(description);
    s += &f!("\n{}\n", voxel_columns());
    s += &"-".repeat(92);
    voxels
        .iter()
        .enumerate()
        .for_each(|(i, v)| s += &f!("\n{i:^6}{}", voxel_row(mesh, v)));
    s
}

#[doc(hidden)]
fn target_heading() -> String {
    f!("{:^72}", "Points to search")
}

#[doc(hidden)]
fn target_columns() -> String {
    let mut s = f!("{:^6}", "id");
    s += &f!("{:^13}", "energy");
    s += &f!("{:^13}", "time");
    s += &f!("{:^13}", "i_coord");
    s += &f!("{:^13}", "j_coord");
    s += &f!("{:^13}", "k_coord");
    s += &f!("{:^6}", "system");
    s
}

#[doc(hidden)]
fn target_row(point: &Point) -> String {
    let mut s = f!("{:^13}", f!("{}", point.e));
    s += &f!("{:^13}", f!("{}", point.t));
    s += &f!("{:^13}", point.i.sci(5, 2));
    s += &f!("{:^13}", point.j.sci(5, 2));
    s += &f!("{:^13}", point.k.sci(5, 2));
    s += &f!(" {:^7}", f!("{}", point.coordinate_type));
    s
}

#[doc(hidden)]
fn voxel_heading(description: &str) -> String {
    let s = f!("Voxels found ({description})");
    f!("{:^92}", s)
}

#[doc(hidden)]
fn voxel_columns() -> String {
    let mut s = f!("{:^6}", "id");
    s += &f!("{:^13}", "energy");
    s += &f!("{:^13}", "time");
    s += &f!("{:^13}", "i_coord");
    s += &f!("{:^13}", "j_coord");
    s += &f!("{:^13}", "k_coord");
    s += &f!("{:^13}", "result");
    s += &f!("{:^8}", "error");
    s
}

#[doc(hidden)]
fn voxel_row(mesh: &Mesh, voxel: &Option<Voxel>) -> String {
    match voxel {
        None => f!("{:^92}", "Not found in mesh"),
        Some(v) => match mesh.voxel_coordinates(v.index) {
            Result::Ok(c) => {
                let mut s = f!("{:^13}", f!("{}", c.energy));
                s += &f!("{:^13}", f!("{}", c.time));
                s += &f!("{:^13}", c.i.sci(5, 2));
                s += &f!("{:^13}", c.j.sci(5, 2));
                s += &f!("{:^13}", c.k.sci(5, 2));
                s += &f!("{:^13}", v.result.sci(5, 2));
                s += &f!("{:^8.4}", v.error);
                s
            }
            Result::Err(_) => {
                let mut s = f!("{:^63}", "Unable to infer voxel coordiantes");
                s += &f!("{:>13}", v.result.sci(5, 2));
                s += &f!("{:>8.4}", v.error);
                s
            }
        },
    }
}
