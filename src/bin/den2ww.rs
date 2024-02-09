// turn off annoying warnings for the example
#![allow(unused_variables)]

// import my library (aka "crate" in Rust)
use meshtal::mesh::{Mesh, Voxel};
use meshtal::read_meshtal_target;
use meshtal::vtk::{mesh_to_vtk, write_vtk, VtkFormat};
use meshtal::weights::mesh_to_ww;
use std::iter::zip;
use itertools::izip;

fn main() {
    println!("Hello, Tim!");

    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ! Read in the three meshes from whatever file format
    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // These will give you a "Mesh" structure with everything ordered
    // to match the column output order, regardless of output format.
    //
    // read_meshtal_target() gives you a "Result", which is either "Ok" if the
    // parsing was a success, or "Err" if there was a problem.
    //
    // The .expect() just gives you the "Mesh" if successful, or exits early
    // with the message string if it failed.
    //
    // variables defined with "let VARIABLE = ", and if you want to be able to
    // change them later then use "let mut VARIABLE = " for mutable.

    // Read the void run mesh tally
    let void_mesh =
        read_meshtal_target("/home/teade/testing/dens2ww/01_void/void.msht", 4).expect("Error: unable to read void mesh");

    // print a summary
    println!("The void mesh summary is:\n{void_mesh}");

    // Read the reduced mesh tally
    let reduced_density_mesh = read_meshtal_target("/home/teade/testing/dens2ww/02_0.1Den_Tot/rd_flux.msht", 4)
        .expect("Error: unable to read reduced_density mesh");

    // Read the uncollided flux mesh
    let uncollided_mesh =
        read_meshtal_target("/home/teade/testing/dens2ww/03_0.1Den_uc_flux/rd_uc_flux.msht", 4).expect("Error: unable to read void mesh");

    // Set gamma parameter
    let gamma = 0.8;
    // Set the density ratio
    let den_ratio = 10.0;


    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ! Now just loop over the voxels and do whatever you want with them
    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // All the voxels are a vector of Voxel objects with their index, result,
    // and error e.g. print out the first 10 voxels
    println!("First 10 voxels in the void mesh:");
    for voxel in &void_mesh.voxels[0..10] {
        println!("{voxel}");
    }

    // Get the build up factors across the mesh
    let buildup = get_buildup(&uncollided_mesh, &reduced_density_mesh, gamma, den_ratio);

    // Get estimate of the forward flux
    let fw_flux = get_fwflux(&buildup, &uncollided_mesh, &void_mesh, den_ratio);

    // Write the foward flux to a vtk file
    let vtk = mesh_to_vtk(&fw_flux);
    write_vtk(vtk.unwrap(), "output.vtk", VtkFormat::Xml).unwrap();

    // Convert the forward flux to a weight window
    let weight_window = mesh_to_ww(&fw_flux, 1.0, 1.0, true);
    weight_window.write("wwinp");
    //

    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ! Example: just add them all together, which I should overload "+" for
    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //let new_results = do_something(&void_mesh, &reduced_density_mesh, &uncollided_mesh);

    println!("\nFirst 10 build up voxels:");
    for voxel in &buildup[0..10] {
        println!("{voxel}");
    }

    println!("\nFirst 10 fw flux voxels:");
    for voxel in &fw_flux.voxels[0..10] {
        println!("{voxel}");
    }

    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ! Then we can chuck this back into a mesh and pass it to the weights
    // ! module for generating wwinp files, or the vtk module for plotting
    // ! the weights/combined mesh results, whatever you want really
    // ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

// Function to work out the buildup factor from the uncollided and total fluxes at reduced density
fn get_buildup(uc_mesh: &Mesh, total_mesh: &Mesh, gamma: f64, den_ratio: f64) -> Vec<Voxel> {

    // start with a copy of the total flux
    let mut build_up = total_mesh.voxels.clone();

    // Loop over voxels to work out buildup factor
    for (bu_voxel, uc_voxel) in zip(&mut build_up, &uc_mesh.voxels) {
        bu_voxel.result = f64::powf(bu_voxel.result/uc_voxel.result,gamma*den_ratio);
    }


    // return the buildup factors for each voxel (note "return" not needed)
    return build_up;

}

// Function to work out an estimate of the foward flux at the full density 
fn get_fwflux(build_up: &Vec<Voxel>, uc_mesh: &Mesh, void_mesh: &Mesh, den_ratio: f64) -> Mesh {

    let mut fw_flux = (*void_mesh).clone();

    // Loop over voxels to work out forward flux
    for (fw_voxel, bu_voxel, uc_voxel, void_voxel) in izip!(&mut fw_flux.voxels, build_up, &uc_mesh.voxels, &void_mesh.voxels) {
        fw_voxel.result = bu_voxel.result * void_voxel.result * f64::powf(uc_voxel.result/void_voxel.result,den_ratio);
    }

    // return the buildup factors for each voxel (note "return" not needed)
    return fw_flux;

}

/// Some function for combining the results and messing with the weights
fn do_something(void: &Mesh, reduced: &Mesh, uncollided: &Mesh) -> Vec<Voxel> {
    // I dunno, for every voxel do (void + reduced) * uncollided
    // details not important, ignoring the error

    // start with a copy of the void voxels
    let mut new_voxels = void.voxels.clone();

    // add the void + reduced voxels manually, damn I need to sort out the + operator
    for (new, r) in zip(&mut new_voxels, &reduced.voxels) {
        new.result += r.result;
    }

    // multiply by uncollided results manually, again sort out a mult operator
    for (new, u) in zip(&mut new_voxels, &uncollided.voxels) {
        new.result *= u.result;
    }

    // return the new voxels (note "return" not needed)
    return new_voxels;
}
