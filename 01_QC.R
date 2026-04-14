library(data.table)
library(EasyQC)
library(stringr)

# ==============================================================================
# Setup paths
# ==============================================================================
input_dir <- "./data/raw_gwas"
output_dir <- "./data/cleaned_gwas"
ref_path_abs <- "./refdata/easyqc/rsmid_machsvs_mapb37.1000G_p3v5.merged_mach_impute.v3.corrpos.gz"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Define EasyQC configuration generation function
# ==============================================================================
create_ecf_content <- function(input_file, short_name, outdir, ref_path) {
  ecf <- paste0(
    "DEFINE	--strSeparator TAB\n",
    "	--strMissing NA\n",
    "	--acolIn CHR;POS;SNP;A1;A2;Effect;SE;Pval;N;INFO\n",
    "	--acolInClasses integer;integer;character;character;character;numeric;numeric;numeric;numeric;numeric\n",
    "	--pathOut ", outdir, "/\n",
    "\n",
    "EASYIN	--fileIn ", input_file, "\n",
    "	--fileInShortName ", short_name, "\n",
    "\n",
    "START EASYQC\n",
    
    # 1. Removal of SNPs with missing values in key fields
    "CLEAN	--rcdClean is.na(CHR) --strCleanName numDrop_Missing_CHR\n",
    "CLEAN	--rcdClean is.na(POS) --strCleanName numDrop_Missing_POS\n",
    "CLEAN	--rcdClean is.na(SNP) --strCleanName numDrop_Missing_SNP\n",
    "CLEAN	--rcdClean is.na(A1) --strCleanName numDrop_Missing_A1\n",
    "CLEAN	--rcdClean is.na(A2) --strCleanName numDrop_Missing_A2\n",
    "CLEAN	--rcdClean is.na(Effect) --strCleanName numDrop_Missing_BETA\n",
    "CLEAN	--rcdClean is.na(SE) --strCleanName numDrop_Missing_SE\n",
    "CLEAN	--rcdClean is.na(Pval) --strCleanName numDrop_Missing_P\n",
    "CLEAN	--rcdClean is.na(N) --strCleanName numDrop_Missing_N\n",
    
    # 2. Removal of SNPs with invalid values
    "CLEAN	--rcdClean (Pval<0)|(Pval>1) --strCleanName numDrop_invalid_P\n",
    "CLEAN	--rcdClean (SE<=0)|(SE==Inf) --strCleanName numDrop_invalid_SE\n",
    "CLEAN	--rcdClean abs(Effect)==Inf --strCleanName numDrop_invalid_BETA\n",
    
    # 3. Exclusion of variants on sex chromosomes (retaining only autosomes)
    "CLEAN	--rcdClean !CHR%in%c(1:22,NA) --strCleanName numDropSNP_ChrXY\n",
    
    # 4. Harmonize Alleles directionality based on the reference (or prep for it)
    "HARMONIZEALLELES --colInA1 A1 --colInA2 A2\n",
    
    # 5. Removal of insertion-deletion variants and multiallelic variants
    "CLEAN	--rcdClean (A1%in%c('I','D')) | (A2%in%c('I','D')) --strCleanName numDrop_INDEL\n",
    "CLEAN	--rcdClean !(length(A1)==1)|!(length(A2)==1) --strCleanName numDrop_multiallelic\n",
    
    # 6. Harmonization of rsIDs against the European 1000 Genomes Phase 3 reference panel
    "CREATECPTID	--fileMap ", ref_path, "\n",
    "		--colMapMarker rsmid\n",
    "		--colMapChr chr\n",
    "		--colMapPos pos\n",
    "		--colInMarker SNP\n",
    "		--colInA1 A1\n",
    "		--colInA2 A2\n",
    "		--colInChr CHR\n",
    "		--colInPos POS\n",
    
    # 7. Removal of duplicate variants based on the generated reference ID (cptid)
    "CLEANDUPLICATES	--colInMarker cptid --strMode removeall\n",
    
    # 8. Add unique ID and finalize output columns
    "ADDCOL	--rcdAddCol paste(CHR, POS, A1, A2, sep=':') --colOut uniqID\n",
    "GETCOLS	--acolOut CHR;POS;SNP;uniqID;cptid;A1;A2;Effect;SE;Pval;N\n",
    
    # 9. Write cleaned summary statistics
    "WRITE	--strPrefix CLEANED_ \n",
    "	--strMissing NA\n",
    "	--strMode txt\n",
    "STOP EASYQC\n"
  )
  return(ecf)
}

# ==============================================================================
# Process files
# ==============================================================================
# You can define your list of input files here (e.g., list.files(input_dir, pattern = "*.txt.gz", full.names = TRUE))
# For example: files <- c("/path/to/migraine_gwas.txt.gz", "/path/to/wm_phenotype1.txt.gz")
files <- list.files(input_dir, pattern = "*.txt.gz", full.names = TRUE)

for (f in files) {
  message("Processing: ", f)
  base <- basename(f)
  dt <- fread(f)
  
  # Map column names dynamically to EasyQC standards
  nms <- names(dt)
  if ("chromosome" %in% nms && !"CHR" %in% nms) setnames(dt, "chromosome", "CHR")
  nms <- names(dt)
  if ("position" %in% nms && !"POS" %in% nms) setnames(dt, "position", "POS")
  nms <- names(dt)
  if ("rs_number" %in% nms && !"SNP" %in% nms) setnames(dt, "rs_number", "SNP")
  nms <- names(dt)
  if (!"SNP" %in% nms && "rsid_ukbb" %in% nms) setnames(dt, "rsid_ukbb", "SNP")
  nms <- names(dt)
  if (!"SNP" %in% nms && "marker" %in% nms) setnames(dt, "marker", "SNP")
  nms <- names(dt)
  if ("other_allele" %in% nms && !"A1" %in% nms) setnames(dt, "other_allele", "A1")
  nms <- names(dt)
  if ("reference_allele" %in% nms && !"A2" %in% nms) setnames(dt, "reference_allele", "A2")
  nms <- names(dt)
  if ("beta" %in% nms && !"Effect" %in% nms) setnames(dt, "beta", "Effect")
  nms <- names(dt)
  if ("se" %in% nms && !"SE" %in% nms) setnames(dt, "se", "SE")
  nms <- names(dt)
  if ("p.value" %in% nms && !"Pval" %in% nms) setnames(dt, "p.value", "Pval")
  
  # Handle -log10 p-values if present
  if (!"Pval" %in% names(dt) && "_-log10_p-value" %in% names(dt)) {
    dt[, Pval := 10^(-as.numeric(get("_-log10_p-value")))]
  }
  
  # Map Sample Size
  if ("n_samples" %in% names(dt) && !"N" %in% names(dt)) {
    dt[, N := as.numeric(n_samples)]
  } else if ("Neff" %in% names(dt) && !"N" %in% names(dt)) {
    dt[, N := as.numeric(Neff)]
  }
  
  # Set dummy INFO if not available
  dt[, INFO := 1]
  
  # Subset to needed columns
  needed <- c("CHR", "POS", "SNP", "A1", "A2", "Effect", "SE", "Pval", "N", "INFO")
  # Keep only columns that exist (in case N is missing in some files)
  needed_present <- intersect(needed, names(dt))
  dt_out <- dt[, ..needed_present]
  
  # Write temporary preQC file
  tmp_file <- file.path(tempdir(), paste0(sub("\\.(tsv|txt)(\\.gz)?$", "", base, ignore.case = TRUE), "_preQC.txt"))
  fwrite(dt_out, tmp_file, sep = "\t", na = "NA", quote = FALSE)
  
  # Create ECF configuration content
  short_name <- sub("\\.(tsv|txt)(\\.gz)?$", "", base, ignore.case = TRUE)
  ecf_content <- create_ecf_content(tmp_file, short_name, output_dir, ref_path_abs)
  
  # Write ECF file
  ecf_file <- file.path(tempdir(), paste0(short_name, ".ecf"))
  writeLines(ecf_content, ecf_file)
  
  # Run EasyQC
  message("Running EasyQC for ", short_name)
  EasyQC(ecf_file)
  
  # Clean up temporary files
  unlink(tmp_file)
  unlink(ecf_file)
}
