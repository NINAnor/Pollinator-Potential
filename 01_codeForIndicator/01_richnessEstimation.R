################################################################################
# Calculate interaction-supported insect indicator
#
# This function combines:
#
#   1. Insect prediction rasters
#   2. Plant prediction rasters
#   3. A plant–insect interaction matrix
#
# to calculate an interaction-supported probability for each insect.
#
# The basic idea is:
#
#   An insect's predicted probability is weighted by the availability/
#   suitability of the plants with which it interacts.
#
# For each insect i:
#
#     Plant support_i =
#
#         sum(Plant probability_j × Interaction strength_ji)
#         ------------------------------------------------
#                sum(Interaction strength_ji)
#
# The final insect indicator is then:
#
#     Interaction-supported insect probability_i =
#
#         Insect probability_i × Plant support_i
#
# The function also handles insect genera that are not present in the
# interaction matrix by assigning them a small interaction weight (`eps`).
#
################################################################################


calculate_interaction_indicator <- function(
    insects,
    plants,
    interactionMatrix,
    eps = 0.001,
    output_dir = "C:/terra_tmp"
) {
  
  
  # ============================================================================
  # 1. Extract insect genus from raster layer names
  # ============================================================================
  #
  # The insect raster layers are assumed to be named using a convention such
  # as:
  #
  #     Genus_species
  #
  # For example:
  #
  #     Bombus_terrestris
  #
  # The code extracts everything before the first underscore:
  #
  #     Bombus_terrestris -> Bombus
  #
  # This is necessary because the interaction matrix is structured at the
  # insect-GENUS level rather than the species level.
  # ============================================================================
  
  insect_genus <- sub(
    "_.*$",
    "",
    names(insects)
  )
  
  
  # ============================================================================
  # 2. Build the plant × insect-genus interaction matrix
  # ============================================================================
  #
  # The interaction matrix is expected to have:
  #
  #     rows    = plant species
  #     columns = insect genera
  #
  # Additional columns such as "...1" and "species" are removed before
  # converting the data frame into a numerical matrix.
  # ============================================================================
  
  M <- as.matrix(
    interactionMatrix[
      ,
      !names(interactionMatrix) %in% c(
        "...1",
        "species"
      ),
      drop = FALSE
    ]
  )
  
  
  # ----------------------------------------------------------------------------
  # Use the plant species names as row names
  # ----------------------------------------------------------------------------
  
  rownames(M) <- interactionMatrix$species
  
  
  # ============================================================================
  # 3. Match plant raster names to interaction-matrix plant names
  # ============================================================================
  #
  # Plant raster names may use underscores:
  #
  #     Achillea_millefolium
  #
  # whereas the interaction matrix may use spaces:
  #
  #     Achillea millefolium
  #
  # Convert underscores to spaces so that the two datasets can be matched.
  # ============================================================================
  
  plant_names <- gsub(
    "_",
    " ",
    names(plants)
  )
  
  
  # ----------------------------------------------------------------------------
  # Check whether every plant prediction raster has a corresponding row in
  # the interaction matrix.
  # ----------------------------------------------------------------------------
  
  missing_plants <- setdiff(
    plant_names,
    rownames(M)
  )
  
  
  # ----------------------------------------------------------------------------
  # Stop the function if any plant is missing.
  #
  # This is an important validation step because otherwise the matrix
  # multiplication below could silently produce incorrect results.
  # ----------------------------------------------------------------------------
  
  if (length(missing_plants) > 0) {
    
    stop(
      "The following plant raster layers are missing ",
      "from the interaction matrix: ",
      paste(
        missing_plants,
        collapse = ", "
      )
    )
  }
  
  
  # ----------------------------------------------------------------------------
  # Reorder the interaction matrix so that its plant rows are in EXACTLY the
  # same order as the plant raster layers.
  #
  # This is critical:
  #
  #     plant raster layer 1 -> interaction matrix row 1
  #     plant raster layer 2 -> interaction matrix row 2
  #     etc.
  # ----------------------------------------------------------------------------
  
  M_plant <- M[
    plant_names,
    ,
    drop = FALSE
  ]
  
  
  # ============================================================================
  # 4. Match insect species to insect genera in the interaction matrix
  # ============================================================================
  #
  # The insect rasters are at species level, whereas the interaction matrix
  # is at genus level.
  #
  # Example:
  #
  #     insect raster:
  #         Bombus_terrestris
  #
  #     extracted genus:
  #         Bombus
  #
  #     interaction matrix:
  #         Bombus
  #
  # `match()` returns the column number in the interaction matrix corresponding
  # to each insect raster layer.
  # ============================================================================
  
  M_index <- match(
    insect_genus,
    colnames(M_plant)
  )
  
  
  # ----------------------------------------------------------------------------
  # Identify which insect genera were successfully matched.
  #
  # TRUE  = genus exists in interaction matrix
  # FALSE = genus is absent from interaction matrix
  # ----------------------------------------------------------------------------
  
  matched <- !is.na(M_index)
  
  
  # ============================================================================
  # 5. Expand the interaction matrix from genus level to species level
  # ============================================================================
  #
  # The original matrix has:
  #
  #     rows    = plant species
  #     columns = insect genera
  #
  # We need:
  #
  #     rows    = plant species
  #     columns = insect species
  #
  # because the final calculation will be performed for every insect raster
  # layer.
  #
  # Every interaction initially receives `eps`.
  #
  # For insects whose genus is present in the interaction matrix, the actual
  # interaction values replace the eps values.
  #
  # For insects whose genus is NOT present, the interaction weights remain eps.
  # ============================================================================
  
  M_insects <- matrix(
    eps,
    nrow = nrow(M_plant),
    ncol = terra::nlyr(insects),
    dimnames = list(
      rownames(M_plant),
      names(insects)
    )
  )
  
  
  # ----------------------------------------------------------------------------
  # Replace the default eps values with actual interaction values for matched
  # insect genera.
  # ----------------------------------------------------------------------------
  
  if (any(matched)) {
    
    M_insects[, matched] <-
      M_plant[
        ,
        M_index[matched],
        drop = FALSE
      ]
  }
  
  
  # ============================================================================
  # 6. Calculate total interaction weight for each insect
  # ============================================================================
  #
  # For each insect:
  #
  #     interaction_weights =
  #         sum of interaction values across all plant species
  #
  # These values are used to normalise the weighted plant support calculation.
  #
  # An insect that interacts with many plants will therefore have a larger
  # total interaction weight than an insect interacting with fewer plants.
  # ============================================================================
  
  interaction_weights <- colSums(
    M_insects
  )
  
  
  # ============================================================================
  # 7. Define a unique temporary output file
  # ============================================================================
  #
  # terra can create temporary files during raster calculations. A unique
  # filename helps prevent file-locking conflicts when the function is run
  # repeatedly.
  #
  # `Sys.getpid()` identifies the current R process.
  #
  # The timestamp provides an additional level of uniqueness.
  # ============================================================================
  
  output_file <- file.path(
    output_dir,
    paste0(
      "plant_support_insects_",
      Sys.getpid(),
      "_",
      format(
        Sys.time(),
        "%Y%m%d_%H%M%S"
      ),
      ".tif"
    )
  )
  
  
  # ----------------------------------------------------------------------------
  # Remove an existing file with the same name, if necessary.
  # ----------------------------------------------------------------------------
  
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  
  
  # ============================================================================
  # 8. Calculate plant support for every insect
  # ============================================================================
  #
  # `terra::lapp()` applies the calculation to each group of raster cells.
  #
  # For every raster cell:
  #
  #     plant_values
  #
  # contains the predicted probability of every plant species.
  #
  # These values are then multiplied by the interaction matrix:
  #
  #     plant_values %*% M_insects
  #
  # producing an interaction-weighted plant support value for every insect.
  #
  # The result is then divided by the total interaction weight of each insect,
  # producing a weighted mean plant probability.
  # ============================================================================
  
  plant_support_insects <- terra::lapp(
    
    plants,
    
    fun = function(...) {
      
      
      # ------------------------------------------------------------------------
      # Extract plant predictions for the current group of cells.
      #
      # Each column corresponds to a plant species.
      # Each row corresponds to a raster cell.
      # ------------------------------------------------------------------------
      
      plant_values <- cbind(...)
      
      
      # ------------------------------------------------------------------------
      # Calculate the interaction-weighted plant support.
      #
      # Matrix multiplication combines:
      #
      #     plant probability
      #
      # with:
      #
      #     plant × insect interaction strength
      #
      # resulting in one value for every insect.
      # ------------------------------------------------------------------------
      
      weighted_sum <-
        plant_values %*% M_insects
      
      
      # ------------------------------------------------------------------------
      # Convert the weighted sum into a weighted mean.
      #
      # Each insect is divided by its total interaction weight.
      #
      # This gives the mean predicted plant probability among the plants
      # associated with that insect, weighted by interaction strength.
      # ------------------------------------------------------------------------
      
      weighted_mean <-
        sweep(
          weighted_sum,
          2,
          interaction_weights,
          "/"
        )
      
      
      # Return the weighted plant-support values for the current cells.
      weighted_mean
    },
    
    # Use a single processing core.
    #
    # This is often safer for large rasters and avoids multiple processes
    # simultaneously writing temporary raster files.
    cores = 1,
    
    # Write the result directly to disk.
    filename = output_file,
    
    overwrite = TRUE
  )
  
  
  # ============================================================================
  # 9. Name the plant-support layers
  # ============================================================================
  #
  # The resulting layers correspond one-to-one with the insect prediction
  # layers.
  #
  # Therefore, give them the same names as the insect rasters.
  # ============================================================================
  
  names(plant_support_insects) <- names(
    insects
  )
  
  
  # ============================================================================
  # 10. Calculate the final interaction-supported insect probability
  # ============================================================================
  #
  # For every insect and every raster cell:
  #
  #     final indicator =
  #
  #         insect prediction
  #         ×
  #         interaction-supported plant probability
  #
  # Therefore, a location receives a high value only when:
  #
  #     1. the insect itself has a high predicted probability, AND
  #
  #     2. the plants associated with that insect also have high predicted
  #        probabilities.
  #
  # This combines species-distribution information with interaction
  # information.
  # ============================================================================
  
  insect_values <-
    insects * plant_support_insects
  
  
  # ============================================================================
  # 11. Return all relevant results
  # ============================================================================
  #
  # Returning the intermediate objects makes it possible to inspect and
  # validate the calculation rather than only retaining the final indicator.
  # ============================================================================
  
  list(
    
    # Final interaction-supported insect probability rasters
    insect_values = insect_values,
    
    # Interaction-supported plant probability for each insect
    plant_support = plant_support_insects,
    
    # Genus extracted from each insect raster layer
    insect_genus = insect_genus,
    
    # TRUE/FALSE indicating whether each insect genus was found
    # in the interaction matrix
    matched = matched,
    
    # Number of matched insect genera
    n_matched = sum(matched),
    
    # Number of unmatched insect genera
    n_unmatched = sum(!matched),
    
    # Small interaction value assigned to unmatched genera
    eps = eps,
    
    # Expanded plant × insect-species interaction matrix
    interaction_matrix = M_insects
  )
}


################################################################################
# Run the interaction indicator calculation
################################################################################
#
# Here:
#
#     eps = 0
#
# means that an insect genus not represented in the interaction matrix receives
# ZERO interaction weight with every plant.
#
# Consequently, unmatched insect genera will have an interaction weight of
# zero, which needs special attention because the weighted-mean calculation
# divides by `interaction_weights`.
################################################################################

weighted_insects_probability <-
  calculate_interaction_indicator(
    insects = insects,
    plants = plants,
    interactionMatrix = interractionMatrix,
    eps = 0.000001
  )
