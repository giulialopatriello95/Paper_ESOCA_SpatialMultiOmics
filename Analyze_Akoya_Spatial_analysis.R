# ==============================================================================

# CosMx spatial transcriptomics analysis

# ==============================================================================

#

# This script performs downstream analysis of a merged CosMx/Seurat object.

#

# Main steps:

# 1. Load required R packages and the processed Seurat object

# 2. Annotate broad cell types

# 3. Subtype B cells according to IgG/IgM expression

# 4. Retrieve and merge spatial coordinates from all fields of view (FOVs)

# 5. Assign individual cells to patients based on sample and spatial coordinates

# 6. Define response groups (CR/NR)

# 7. Calculate spatial cell-type neighbourhoods using k-nearest neighbours

# 8. Summarize cell-type contacts for each patient

#

# Input:

# - complete.Rdata

# Contains the processed merged Seurat object named `merged_seurat`.

#

# Output:

# - all_patients_contact_with_other_bcell_k20.csv

# Patient-level summary of cell-type neighbourhood composition.

#

# ==============================================================================

# ------------------------------------------------------------------------------

# 1. Load required packages

# ------------------------------------------------------------------------------

library(Seurat)
library(InSituCor)
library(dplyr)
library(reshape)
library(mclust)
library(ggplot2)
library(reshape2)
library(magrittr)
library(future)
library(RANN)

# Configure the future framework.

# Sequential execution is used to avoid excessive memory usage when working

# with large spatial transcriptomics objects.

options(future.globals.maxSize = 1e9)
plan("sequential")

# ------------------------------------------------------------------------------

# 2. Load processed Seurat object

# ------------------------------------------------------------------------------

# Load the previously processed and integrated Seurat object.

load("complete.Rdata")

# ------------------------------------------------------------------------------

# 3. Assign broad cell-type annotations

# ------------------------------------------------------------------------------

# Initialize all cells as "Other".

merged_seurat$cell_type <- "Other"

# Assign broad cell types according to the existing Seurat cluster identities.

merged_seurat$cell_type[
Idents(merged_seurat) %in% c("3", "4", "6", "11", "13")
] <- "Tumor"

merged_seurat$cell_type[
Idents(merged_seurat) %in% c("10")
] <- "Dentritic cells"

merged_seurat$cell_type[
Idents(merged_seurat) %in% c("2", "8")
] <- "B cells"

merged_seurat$cell_type[
Idents(merged_seurat) %in% c("0", "7")
] <- "T cells"

# ------------------------------------------------------------------------------

# 4. Identify IgG+ and IgM+ B-cell populations

# ------------------------------------------------------------------------------

# Extract B cells from the complete dataset.

B_cells <- subset(
merged_seurat,
subset = cell_type == "B cells"
)

# Use IgM and IgG expression measurements as variable features for

# visualization and downstream clustering of the B-cell population.

VariableFeatures(B_cells) <- c(
"Cell..IgM.mean",
"Cell..IgG.mean"
)

# Construct a nearest-neighbour graph using the integrated RPCA reduction.

B_cells <- FindNeighbors(
B_cells,
reduction = "integrated.rpca",
dims = 1:3
)

# Perform clustering of the B-cell population.

B_cells <- FindClusters(
B_cells,
resolution = 0.6
)

# Generate a UMAP representation of B-cell subclusters.

B_cells <- RunUMAP(
B_cells,
dims = 1:3
)

DimPlot(
B_cells,
label = TRUE
)

# Identify cells belonging to clusters characterized by IgG expression.

IgG_cells <- rownames(
[B_cells@meta.data](mailto:B_cells@meta.data)[
Idents(B_cells) %in% c("8", "4", "10"),
]
)

# Identify cells belonging to clusters characterized by IgM expression.

IgM_cells <- rownames(
[B_cells@meta.data](mailto:B_cells@meta.data)[
Idents(B_cells) %in% c("3", "7"),
]
)

# Add the B-cell subtypes back to the complete Seurat object.

merged_seurat$cell_type[
rownames([merged_seurat@meta.data](mailto:merged_seurat@meta.data)) %in% IgG_cells
] <- "IgG+ B cells"

merged_seurat$cell_type[
rownames([merged_seurat@meta.data](mailto:merged_seurat@meta.data)) %in% IgM_cells
] <- "IgM+ B cells"

# ------------------------------------------------------------------------------

# 5. Retrieve spatial coordinates from all FOVs

# ------------------------------------------------------------------------------

# Extract cell coordinates from each CosMx field of view.

coords_1LB  <- GetTissueCoordinates(merged_seurat, image = "X1LB")
coords_3LB  <- GetTissueCoordinates(merged_seurat, image = "X3LB")
coords_6LB  <- GetTissueCoordinates(merged_seurat, image = "X6LB")
coords_7LB  <- GetTissueCoordinates(merged_seurat, image = "X7LB")
coords_9LB  <- GetTissueCoordinates(merged_seurat, image = "X9LB")
coords_11LB <- GetTissueCoordinates(merged_seurat, image = "X11LB")
coords_13LB <- GetTissueCoordinates(merged_seurat, image = "X13LB")
coords_15LB <- GetTissueCoordinates(merged_seurat, image = "X15LB")
coords_17LB <- GetTissueCoordinates(merged_seurat, image = "X17LB")
coords_19LB <- GetTissueCoordinates(merged_seurat, image = "X19LB")
coords_21LB <- GetTissueCoordinates(merged_seurat, image = "X21LB")
coords_23LB <- GetTissueCoordinates(merged_seurat, image = "X23LB")
coords_25LB <- GetTissueCoordinates(merged_seurat, image = "X25LB")
coords_27LB <- GetTissueCoordinates(merged_seurat, image = "X27LB")

# Coordinates from additional CosMx samples/FOVs.

coords_7599  <- GetTissueCoordinates(
merged_seurat,
image = "X7599_part2"
)

coords_75293 <- GetTissueCoordinates(
merged_seurat,
image = "X00075293"
)

coords_75298 <- GetTissueCoordinates(
merged_seurat,
image = "X00075298"
)

coords_75292 <- GetTissueCoordinates(
merged_seurat,
image = "X00075292"
)

# Combine coordinates from all FOVs into a single data frame.

all_coords <- do.call(
rbind,
list(
coords_1LB,
coords_3LB,
coords_6LB,
coords_7LB,
coords_9LB,
coords_11LB,
coords_13LB,
coords_15LB,
coords_17LB,
coords_19LB,
coords_21LB,
coords_23LB,
coords_25LB,
coords_27LB,
coords_7599,
coords_75293,
coords_75298,
coords_75292
)
)

# ------------------------------------------------------------------------------

# 6. Merge spatial coordinates with cell metadata

# ------------------------------------------------------------------------------

# Extract Seurat metadata and add cell IDs as an explicit column.

metadata <- [merged_seurat@meta.data](mailto:merged_seurat@meta.data)
metadata$cell <- rownames(metadata)

# Merge cell metadata with the spatial coordinates using the cell ID.

metadata_merged <- merge(
metadata,
all_coords,
by.x = "cell",
by.y = "cell",
all.x = TRUE,
no.dups = TRUE
)

# ------------------------------------------------------------------------------

# 7. Assign patient IDs based on sample and spatial coordinates

# ------------------------------------------------------------------------------

# Initialize the patient column.

metadata_merged$Patient <- ""

# Patient assignment is based on the CosMx sample/FOV and the spatial

# coordinate boundaries separating the individual patient regions.

metadata_merged$Patient[
metadata_merged$sample == "1LB" &
metadata_merged$y < 10000
] <- "CR8459"

metadata_merged$Patient[
metadata_merged$sample == "1LB" &
metadata_merged$y > 10000
] <- "NR3014"

metadata_merged$Patient[
metadata_merged$sample == "3LB" &
metadata_merged$y < 10000
] <- "CR6290"

metadata_merged$Patient[
metadata_merged$sample == "3LB" &
metadata_merged$y > 10000
] <- "NR1231"

metadata_merged$Patient[
metadata_merged$sample == "6LB" &
metadata_merged$y < 10000
] <- "CR6858"

metadata_merged$Patient[
metadata_merged$sample == "6LB" &
metadata_merged$y > 10000
] <- "NR2591"

metadata_merged$Patient[
metadata_merged$sample == "7LB" &
metadata_merged$y < 10000
] <- "CCR3440"

metadata_merged$Patient[
metadata_merged$sample == "7LB" &
metadata_merged$y > 10000
] <- "NR3497"

metadata_merged$Patient[
metadata_merged$sample == "9LB" &
metadata_merged$y < 10000
] <- "CR9404"

metadata_merged$Patient[
metadata_merged$sample == "9LB" &
metadata_merged$y > 10000
] <- "NR0744"

metadata_merged$Patient[
metadata_merged$sample == "11LB" &
metadata_merged$y < 10000
] <- "CR1764"

metadata_merged$Patient[
metadata_merged$sample == "11LB" &
metadata_merged$y > 10000
] <- "NR4406"

metadata_merged$Patient[
metadata_merged$sample == "13LB" &
metadata_merged$y < 10000
] <- "CR6004"

metadata_merged$Patient[
metadata_merged$sample == "13LB" &
metadata_merged$y > 10000
] <- "NR7418"

metadata_merged$Patient[
metadata_merged$sample == "15LB" &
metadata_merged$y < 10000
] <- "CR6563"

metadata_merged$Patient[
metadata_merged$sample == "15LB" &
metadata_merged$y > 10000
] <- "NR5037"

metadata_merged$Patient[
metadata_merged$sample == "17LB" &
metadata_merged$y < 10000
] <- "CR1034"

metadata_merged$Patient[
metadata_merged$sample == "17LB" &
metadata_merged$y > 10000
] <- "NR6466"

metadata_merged$Patient[
metadata_merged$sample == "19LB" &
metadata_merged$y < 10000
] <- "CR2410"

metadata_merged$Patient[
metadata_merged$sample == "19LB" &
metadata_merged$y > 10000
] <- "NR2875"

metadata_merged$Patient[
metadata_merged$sample == "21LB" &
metadata_merged$y < 10000
] <- "CR9418"

metadata_merged$Patient[
metadata_merged$sample == "21LB" &
metadata_merged$y > 10000
] <- "NR2059"

metadata_merged$Patient[
metadata_merged$sample == "23LB" &
metadata_merged$y < 10000
] <- "CR4291"

metadata_merged$Patient[
metadata_merged$sample == "23LB" &
metadata_merged$y > 10000
] <- "NR7252"

metadata_merged$Patient[
metadata_merged$sample == "25LB" &
metadata_merged$y < 10000
] <- "CR4423"

metadata_merged$Patient[
metadata_merged$sample == "25LB" &
metadata_merged$y > 10000
] <- "NR6833"

metadata_merged$Patient[
metadata_merged$sample == "27LB" &
metadata_merged$y < 10000
] <- "CR5726"

metadata_merged$Patient[
metadata_merged$sample == "27LB" &
metadata_merged$y > 10000
] <- "NR8347"

# Additional samples containing multiple patient regions.

metadata_merged$Patient[
metadata_merged$sample == "00075293" &
metadata_merged$y > 10000
] <- "NR21-I-11925"

metadata_merged$Patient[
metadata_merged$sample == "00075293" &
metadata_merged$y < 10000
] <- "CR21-I-11612"

metadata_merged$Patient[
metadata_merged$sample == "00075292" &
metadata_merged$y > 15000
] <- "CR52010783"

metadata_merged$Patient[
metadata_merged$sample == "00075292" &
metadata_merged$y < 15000 &
metadata_merged$y > 5000
] <- "NR51929838"

metadata_merged$Patient[
metadata_merged$sample == "00075292" &
metadata_merged$y < 5000
] <- "NR51800155"

metadata_merged$Patient[
metadata_merged$sample == "00075298" &
metadata_merged$y < 5000
] <- "CR51806411"

metadata_merged$Patient[
metadata_merged$sample == "00075298" &
metadata_merged$y > 5000 &
metadata_merged$y < 15000
] <- "NR21-I-28696"

metadata_merged$Patient[
metadata_merged$sample == "00075298" &
metadata_merged$y > 15000
] <- "NR22-I-09272"

metadata_merged$Patient[
metadata_merged$sample == "7599_part2" &
metadata_merged$y < 5000
] <- "CR24-I-00321"

metadata_merged$Patient[
metadata_merged$sample == "7599_part2" &
metadata_merged$y > 5000 &
metadata_merged$y < 15000
] <- "NR24-I-01034"

metadata_merged$Patient[
metadata_merged$sample == "7599_part2" &
metadata_merged$y > 15000
] <- "NR24-I-14933"

# Restore the original cell IDs as row names.

rownames(metadata_merged) <- metadata_merged$cell

# Add patient information to the Seurat object.

merged_seurat <- AddMetaData(
merged_seurat,
metadata = metadata_merged[, "Patient", drop = FALSE]
)

# ------------------------------------------------------------------------------

# 8. Define response group

# ------------------------------------------------------------------------------

# Initialize all samples as non-responders.

merged_seurat$Condition <- "NR"

# Assign CR to patients whose IDs contain "CR".

merged_seurat$Condition[
grepl("CR", [merged_seurat@meta.data](mailto:merged_seurat@meta.data)$Patient)
] <- "CR"

# ------------------------------------------------------------------------------

# 9. Custom Seurat subsetting function

# ------------------------------------------------------------------------------

# This function extends the standard Seurat subsetting workflow to objects

# containing multiple FOVs. It ensures that cells are correctly retained

# across FOV-specific coordinate information.

subset_opt <- function(
object = NULL,
subset,
cells = NULL,
idents = NULL,
Update.slots = TRUE,
Update.object = TRUE,
...
) {

if (Update.slots) {
message("Updating object slots...")
object %<>% UpdateSlots()
}

message("Cloning object...")
obj_subset <- object

# If integer indices are supplied, convert them to cell IDs.

if (all(is.integer(cells))) {
cells <- Cells(obj_subset)[cells]
}

if (!missing(subset) || !is.null(idents)) {
message(
"Extracting cells matched to `subset` and/or `idents`"
)
}

# Handle FOV objects separately.

if (class(obj_subset) == "FOV") {

```
message("Object class is FOV")
cells <- Cells(obj_subset)
```

} else if (
!class(obj_subset) == "FOV" &&
!missing(subset)
) {

```
subset <- enquo(arg = subset)

cells <- WhichCells(
  object = obj_subset,
  cells = cells,
  idents = idents,
  expression = subset,
  return.null = TRUE,
  ...
)
```

} else if (
!class(obj_subset) == "FOV" &&
!is.null(idents)
) {

```
cells <- WhichCells(
  object = obj_subset,
  cells = cells,
  idents = idents,
  return.null = TRUE,
  ...
)
```

} else if (is.null(cells)) {

```
cells <- Cells(obj_subset)
```

}

# Check whether selected cells are present in each FOV.

if (class(obj_subset) == "FOV") {

```
message("Matching cells for FOV object...")
cells_check <- any(
  obj_subset %>% Cells %in% cells
)
```

} else {

```
message("Matching cells in FOVs...")

cells_check <- lapply(
  Images(obj_subset) %>% seq,
  function(i) {
    any(
      obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>%
        Cells %in% cells
    )
  }
) %>%
  unlist()
```

}

# If all selected cells are present in all FOVs, standard subsetting

# can be performed.

if (all(cells_check)) {

```
message(
  "Cell subsets are found in all FOVs!",
  "\nSubsetting object..."
)

obj_subset %<>%
  base::subset(
    cells = cells,
    idents = idents,
    ...
  )
```

} else {

```
# If selected cells are present only in some FOVs, subset each
# relevant FOV individually.
fovs <- lapply(
  Images(obj_subset) %>% seq,
  function(i) {

    if (
      any(
        obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>%
          Cells %in% cells
      )
    ) {

      message(
        "Subsetting FOV: ",
        Images(obj_subset)[i]
      )

      base::subset(
        x = obj_subset[[Images(obj_subset)[i]]],
        cells = cells,
        idents = idents,
        ...
      )
    }
  }
)

# Replace the original FOVs with their subsetted versions.
for (i in fovs %>% seq) {
  obj_subset[[Images(object)[i]]] <- fovs[[i]]
}
```

}

# Perform final cell-level subsetting.

obj_subset %<>%
base::subset(cells = cells, ...)

# Update the Seurat object after subsetting.

if (
Update.object &&
!class(obj_subset) == "FOV"
) {

```
message("Updating object...")
obj_subset %<>%
  UpdateSeuratObject()
```

}

message("Object is ready!")

return(obj_subset)
}

# ------------------------------------------------------------------------------

# 10. Prepare Seurat object for spatial analysis

# ------------------------------------------------------------------------------

# Store the sample/FOV identifier as the FOV metadata field.

merged_seurat$FOV <- merged_seurat$sample

# Use the RNA assay for downstream analysis.

DefaultAssay(merged_seurat) <- "RNA"

# Remove the additional assay if it is no longer required.

merged_seurat[["otherAssay"]] <- NULL

# ------------------------------------------------------------------------------

# 11. Define patient list

# ------------------------------------------------------------------------------

# Extract unique patient IDs and remove cells without a patient assignment.

patient_ids <- unique(
[merged_seurat@meta.data](mailto:merged_seurat@meta.data)$Patient
)

patient_ids <- patient_ids[
patient_ids != ""
]

# Initialize a list to store the results for each patient.

all_patient_results <- list()

# ------------------------------------------------------------------------------

# 12. Spatial cell-cell neighbourhood analysis

# ------------------------------------------------------------------------------

# Number of nearest neighbours considered for each cell.

k <- 20

# Loop through all patients independently.

for (patient_id in patient_ids) {

message(
"Processing patient: ",
patient_id
)

# --------------------------------------------------------------------------

# 12.1 Subset the Seurat object to the current patient

# --------------------------------------------------------------------------

CosMX_patient <- subset_opt(
merged_seurat,
Patient == patient_id
)

# --------------------------------------------------------------------------

# 12.2 Extract spatial coordinates

# --------------------------------------------------------------------------

xy <- GetTissueCoordinates(
CosMX_patient
)

# Shift coordinates so that the minimum x and y coordinates are zero.

xy$x <- xy$x - min(xy$x)
xy$y <- xy$y - min(xy$y)

rownames(xy) <- xy$cell

coords <- xy[, c("x", "y")]
cells <- rownames(coords)

# --------------------------------------------------------------------------

# 12.3 Identify nearest neighbours

# --------------------------------------------------------------------------

# Find the k+1 nearest neighbours for each cell.

# The additional neighbour corresponds to the cell itself and is removed

# below.

nn <- nn2(
coords,
k = k + 1
)

# Convert neighbour indices into cell IDs.

nn_cells <- apply(
nn$nn.idx[, -1],
2,
function(x) cells[x]
)

rownames(nn_cells) <- cells

# --------------------------------------------------------------------------

# 12.4 Generate all cell-type combinations

# --------------------------------------------------------------------------

# Identify the cell types present in the current patient.

celltype_names <- unique(
as.character(
CosMX_patient$cell_type
)
)

# Generate every possible ordered pair of cell types.

feat_pairs <- expand.grid(
A = celltype_names,
B = celltype_names,
stringsAsFactors = FALSE
)

# Container for patient-level results.

m_all <- data.frame(
a = character(),
b = character(),
abs_count = numeric(),
prop = double(),
patient = character()
)

# --------------------------------------------------------------------------

# 12.5 Calculate cell-type neighbourhood composition

# --------------------------------------------------------------------------

for (i in seq_len(nrow(feat_pairs))) {

```
a <- feat_pairs[i, 1]
b <- feat_pairs[i, 2]

# Identify cells belonging to cell type A and B.
a_names <- rownames(
  CosMX_patient@meta.data[
    CosMX_patient$cell_type == a,
    ,
    drop = FALSE
  ]
)

b_names <- rownames(
  CosMX_patient@meta.data[
    CosMX_patient$cell_type == b,
    ,
    drop = FALSE
  ]
)

if (
  length(a_names) > 0 &&
  length(b_names) > 0
) {

  # For each cell of type B, count the number and proportion of
  # neighbouring cells belonging to cell type A.
  res <- t(
    sapply(
      b_names,
      function(cell) {

        neighs <- nn_cells[cell, ]

        abs_count <- sum(
          CosMX_patient@meta.data[
            neighs,
            "cell_type"
          ] == a
        )

        prop <- abs_count / k

        c(
          abs_count = abs_count,
          proportion = prop
        )
      }
    )
  )

  res_df <- as.data.frame(res)

  res_df$a <- a
  res_df$b <- b
  res_df$patient <- patient_id

  # Average the number and proportion of A cells across all B cells
  # within the current patient.
  new_row <- data.frame(
    a = a,
    b = b,
    abs_count = mean(
      res_df$abs_count
    ),
    prop = mean(
      res_df$proportion
    ),
    patient = patient_id
  )

} else {

  # If either cell type is absent, record zero contact.
  new_row <- data.frame(
    a = a,
    b = b,
    abs_count = 0,
    prop = 0,
    patient = patient_id
  )
}

# Append the current cell-type pair to the patient-level results.
m_all <- rbind(
  m_all,
  new_row
)
```

}

# Store results for the current patient.

all_patient_results[[patient_id]] <- m_all

# Remove temporary objects and release memory.

rm(
CosMX_patient,
xy,
coords,
cells,
m_all
)

gc()
}

# ------------------------------------------------------------------------------

# 13. Combine patient-level results

# ------------------------------------------------------------------------------

# Combine all patient-level neighbourhood results into a single data frame.

combined_results <- do.call(
rbind,
all_patient_results
)

# ------------------------------------------------------------------------------

# 14. Save results

# ------------------------------------------------------------------------------

write.csv(
combined_results,
"all_patients_contact_with_other_bcell_k20.csv",
row.names = FALSE
)

# ==============================================================================

# End of script

# ==============================================================================
