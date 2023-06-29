# //////////////////////////////////////////////////////////////////////////////
docs_main = dict(
    main = """
  +===================================+
  | Computational Biophyisics Project |
  +===================================+""",
    calc = "precalculate the values required for the different visualizations",
    bse = "visualize the BSE line plots",
    cluster = "visualize the hierarchical clustering interactive scatter plots",
    cmap0 = "visualize the CMAP heatmap plots",
    cmap1 = "visualize the DCMAP and DCMD interactive heatmap plots",
    cmap2 = "visualize the active site CMAP line plots",
    pca = "visualize the PCA interactive scatter plots, including the cumulative variance of the systems",
    pyinteraph = "visualize the Pyinteraph heatmap plots",
    rama = "visualize the RAMA interactive scatter plots",
    rgyr = "visualize the RGYR scatter plots",
    rmsd0 = "visualize the RMSD heatmap plots. WARNING: this option is memory intensive, be sure to have at least 2.6 GB of RAM available",
    rmsd1 = "visualize the RMSD interactive scatter plots",
    rmsf0 = "visualize the RMSF line plots (comparison between runs of the same system)",
    rmsf1 = "visualize the RMSF line plots (comparison between different systems)",
    sasa = "visualize the SASA scatter plots",
)

# //////////////////////////////////////////////////////////////////////////////
docs_wald = dict(
    main = """
  +===============================================+
  | Water Average Local Densities (WALD) pipeline |
  +===============================================+
    * Go through step 0 by selecting the reference frame of interest and the allowed RMSD treshold of "similarity".
    * Execute steps 1 and 2.
    * BEFORE step 3:
      + Open the VMD file generated in step 1.
      + Reset the view and toggle off the axes.
      + Select the reference frame as indicated in the console.
      + File -> Render... -> "Render the current scene using:" Wavefront (OBJ and MTL).
      + Rename the file with the same run prefix as the VMD file, but with an OBJ extension (e.g. mt2_rep0.obj).
      + Start Rendering.
    * Execute step 3.
    * Visualize the result with step 4 as many times as needed.
""",
    rmsd = "select a set of frames that are structurally similar to a reference frame, defined by having RMSD values lower than a given treshold",
    calc = "calculate the 3D WALD grid",
    vmd = "generate the VMD visualization state required to render the meshes of the protein representation",
    meshes = "parse the protein meshes and convert them to a more appropriate format for Plotly",
    plot = "visualize the protein and WALD meshgrids interactively",
)

# //////////////////////////////////////////////////////////////////////////////
