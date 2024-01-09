//! Module for dealing with UKAEA CuV posvol binaries
//!
//! Very useful information for the dominant cells in each voxels, broken down
//! into what I will be calling `sub-voxels`.
//!
//! Foe example, a resolution of 5x5x5 on the CuV IDUM card will break every
//! voxel up into 125 regions and sample each region to find the dominant cell.
//!
//! This knowledge can be used for much finer resolution plotting in VTK outputs
//! for CuV meshes.
//!
//! ## Reading a posvol file
//!
//! To read the binary file, the path is provided via the readers module.
//!
//! ```rust
//! # use meshtal::readers;
//! let posvol = readers::read_posvol_file("./data/posvol/plot_fmesh_614.bin").unwrap();
//! ```  
//!
//! This deserialises from the binary format into useful data structures.
//!

// internal modules
use crate::utils::f;

// external crates
use serde::{Deserialize, Serialize};

/// Representation og data in a posvol file
///
/// The byte layout is very simple. The 6 dimension values in the first block
/// are stored as [Dimensions].
///
/// ```text
/// <block byte length>
///     <resolution i> <resolution j> <resolution k>
///     <iints+1> <jints+1> <kints+1>
/// <block byte length>
/// ```
///
/// The second block contains all cell data in a continuous array, and is stored
/// as a vector of cell values.
///
/// ```text
/// <block byte length>
///     <voxel 0, subvoxel 0> <voxel 0, subvoxel 1>  <voxel 0, subvoxel 2> ...
///     <voxel 1, subvoxel 0> <voxel 1, subvoxel 1>  <voxel 1, subvoxel 2> ...
///     ... and so on
/// <block byte length>
/// ```
#[derive(Debug, Serialize)]
pub struct Posvol {
    /// The dimensions given in the first block of data
    pub dimensions: Dimensions,
    /// List of dimonant cells for every subvoxel
    pub cells: Vec<i32>,
}

impl Posvol {
    /// Vector of subvoxel cell groups
    ///
    /// Extremely common to iterate over the voxels in chunks of subvoxel cells.
    pub fn subvoxels(&self) -> Vec<&[i32]> {
        self.cells
            .chunks_exact(self.dimensions.number_of_subvoxels())
            .collect()
    }

    /// Number of voxels expected in the file
    pub fn number_of_voxels(&self) -> usize {
        self.dimensions.number_of_voxels()
    }

    /// Number of samples per voxel expected in the file
    pub fn number_of_subvoxels(&self) -> usize {
        self.dimensions.number_of_subvoxels()
    }

    /// Total number of cells expected in the file
    pub fn number_of_cells(&self) -> usize {
        self.dimensions.number_of_cells()
    }
}

impl std::fmt::Display for Posvol {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut s = "Posvol {\n".to_string();
        s += &f!(
            "    voxels: {} ({}x{}x{})\n",
            self.dimensions.number_of_voxels(),
            self.dimensions.n_x - 1,
            self.dimensions.n_y - 1,
            self.dimensions.n_z - 1
        );
        s += &f!(
            "    subvoxels: {} ({}x{}x{})\n",
            self.dimensions.number_of_subvoxels(),
            self.dimensions.res_x,
            self.dimensions.res_y,
            self.dimensions.res_z
        );
        s += &f!(
            "    cells: {} ({}x{})\n}}",
            self.dimensions.number_of_cells(),
            self.dimensions.number_of_voxels(),
            self.dimensions.number_of_subvoxels(),
        );

        write!(f, "{}", s)
    }
}

/// Dimension values in the first [Posvol] data block
///
/// The 6 dimension values in the first block of data are stored here.
///
/// Not intended to be edited but do what you want with it. Fields correspond to
/// the sample resolution in x, y, and z dimensions, and the number of mesh
/// bounds in each mesh coordinate axis.
#[derive(Deserialize, Serialize, Debug)]
pub struct Dimensions {
    pub res_x: i32,
    pub res_y: i32,
    pub res_z: i32,
    pub n_x: i32,
    pub n_y: i32,
    pub n_z: i32,
}

impl Dimensions {
    /// Number of voxels expected in the file
    pub fn number_of_voxels(&self) -> usize {
        ((self.n_x - 1) * (self.n_y - 1) * (self.n_z - 1)) as usize
    }

    /// Number of samples per voxel expected in the file
    pub fn number_of_subvoxels(&self) -> usize {
        (self.res_x * self.res_y * self.res_z) as usize
    }

    /// Total number of cells expected in the file
    pub fn number_of_cells(&self) -> usize {
        self.number_of_voxels() * self.number_of_subvoxels()
    }

    /// Expected size of full cells array
    pub fn cell_array_byte_length(&self) -> i32 {
        (self.number_of_cells() * std::mem::size_of::<i32>()) as i32
    }
}

impl std::fmt::Display for Dimensions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}
