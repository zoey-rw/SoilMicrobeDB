# Download Mycocosm fungal genomes for Soil Microbe Database
# Uses RSelenium to automate browser-based downloads (requires authentication)
#
# Usage:
#   Rscript download_actual_myco_genomes.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_actual_myco_genomes.r  # Test with 10 genomes
#
# Environment variables:
#   MYCOCOSM_EMAIL: JGI Mycocosm login email (required)
#   MYCOCOSM_PASSWORD: JGI Mycocosm login password (required)
#   RSELENIUM_PORT: Port for RSelenium server (default: 5003)
#   RSELENIUM_BROWSER: Browser to use (default: firefox)

library(httr)
library(jsonlite)
library(RSelenium)
library(rvest)
library(tidyverse)

# Load centralized configuration
# Find the create_database directory (go up from process_genomes/Mycocosm)
script_dir <- tryCatch({
  # Try to get script path from commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    dirname(normalizePath(script_path))
  } else {
    # Fallback: try sys.frame
    dirname(normalizePath(sys.frame(1)$ofile))
  }
}, error = function(e) {
  getwd()
})
if (length(script_dir) == 0 || script_dir == "." || is.na(script_dir)) {
  script_dir <- getwd()
}
# Go up to create_database directory (from process_genomes/Mycocosm)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get Mycocosm-specific configuration
mycocosm_config <- get_source_config("Mycocosm")
source_name <- "Mycocosm"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting Mycocosm genome download", source_name = source_name)

# Check for required credentials
mycocosm_email <- Sys.getenv("MYCOCOSM_EMAIL", unset = NA)
mycocosm_password <- Sys.getenv("MYCOCOSM_PASSWORD", unset = NA)

if (is.na(mycocosm_email) || is.na(mycocosm_password)) {
  stop("MYCOCOSM_EMAIL and MYCOCOSM_PASSWORD environment variables are required.\n",
       "Please set them before running this script:\n",
       "  export MYCOCOSM_EMAIL='your_email@example.com'\n",
       "  export MYCOCOSM_PASSWORD='your_password'")
}

# RSelenium configuration
rselenium_port <- as.integer(Sys.getenv("RSELENIUM_PORT", unset = "5003"))
rselenium_browser <- Sys.getenv("RSELENIUM_BROWSER", unset = "firefox")

# Determine output directory (prefer local, fallback to remote)
output_dir <- mycocosm_config$local_genome_dir
if (!dir.exists(output_dir)) {
  # Try remote directory
  if (dir.exists(mycocosm_config$remote_genome_dir)) {
    output_dir <- mycocosm_config$remote_genome_dir
    log_message(paste("Using remote directory:", output_dir), source_name = source_name)
  } else {
    # Create local directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    log_message(paste("Created local directory:", output_dir), source_name = source_name)
  }
} else {
  log_message(paste("Using local directory:", output_dir), source_name = source_name)
}

# Download catalog
log_message("Downloading Mycocosm catalog", source_name = source_name)
tryCatch({
  mycocosm_data_in <- read.csv(
    mycocosm_config$catalog_url,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    # Rename columns to match expected names (handle spaces and special characters)
    rename(
      row = `##`,
      Name = name,
      ID = portal,
      NCBI_TaxID = `NCBI Taxon`,
      gene_count = `#of genes`,
      is_restricted = `is restricted`,
      is_public = `is public`,
      is_published = `is published`,
      is_superseded = `is superseded`,
      `superseded by` = `superseded by`,
      publications = `publication(s)`,
      pubmed_id = `pubmed id(s)`,
      doi_id = `doi id(s)`
    ) %>%
    mutate(NCBI_TaxID = as.numeric(NCBI_TaxID))
  mycocosm_data_in$Name <- gsub('\\"', "", mycocosm_data_in$Name)
  mycocosm_data_in$Go.Download.Link <- str_c(
    'https://genome.jgi.doe.gov/portal/', 
    mycocosm_data_in$ID, 
    "/download/", 
    mycocosm_data_in$ID, 
    "_AssemblyScaffolds_Repeatmasked.fasta.gz"
  )
  log_message(paste("Downloaded catalog with", nrow(mycocosm_data_in), "genomes"), source_name = source_name)
}, error = function(e) {
  stop(paste("Failed to download Mycocosm catalog:", e$message))
})

# Filter to published genomes
download_df <- mycocosm_data_in %>% filter(is_published == "Y")
log_message(paste("Found", nrow(download_df), "published genomes"), source_name = source_name)

# Check for already downloaded genomes
log_message("Checking for already downloaded genomes", source_name = source_name)
list_already_downloaded <- list.files(
  output_dir,
  pattern = mycocosm_config$genome_pattern
)

# Also check remote directory if different
if (output_dir != mycocosm_config$remote_genome_dir && dir.exists(mycocosm_config$remote_genome_dir)) {
  remote_downloaded <- list.files(
    mycocosm_config$remote_genome_dir,
    pattern = mycocosm_config$genome_pattern
  )
  list_already_downloaded <- unique(c(list_already_downloaded, remote_downloaded))
  log_message(paste("Found", length(remote_downloaded), "additional genomes in remote directory"), source_name = source_name)
}

log_message(paste("Found", length(list_already_downloaded), "already downloaded genomes"), source_name = source_name)

# Subset to only those missing from our downloaded files
download_df_subset <- download_df[
  which(!basename(download_df$Go.Download.Link) %in% list_already_downloaded),
]

log_message(paste("Need to download", nrow(download_df_subset), "genomes"), source_name = source_name)

# Apply test mode limit if specified
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes"), source_name = source_name)
  download_df_subset <- download_df_subset %>% slice_head(n = MAX_GENOMES)
}

if (nrow(download_df_subset) == 0) {
  log_message("All genomes already downloaded. Exiting.", source_name = source_name)
  log_message("Mycocosm download script complete", source_name = source_name)
  quit(save = "no", status = 0)
}

# Initialize RSelenium
log_message("Initializing RSelenium", source_name = source_name)
fprof <- makeFirefoxProfile(list(
  browser.download.manager.showWhenStarting = FALSE,
  browser.download.dir = normalizePath(output_dir),
  browser.helperApps.neverAsk.saveToDisk = "text/csv.fasta.gz",
  browser.download.folderList = 2L
))

rd <- rsDriver(
  browser = rselenium_browser,
  chromever = NULL,
  extraCapabilities = fprof,
  port = rselenium_port
)
remDr <- rd$client

# Function to download genomes using RSelenium
download_from_links <- function(run.list, save.folder) {
  run.length <- nrow(run.list)
  
  if (run.length == 0) {
    return()
  }
  
  # Ensure browser closes even if function errors
  on.exit({
    tryCatch({
      remDr$close()
    }, error = function(e) {
      # Browser may already be closed, ignore
    })
  }, add = TRUE)
  
  # Start the remote driver
  remDr$open()
  remDr$screenshot(TRUE)
  
  # Navigate to initial link to prompt login page
  initial.link <- 'https://genome.jgi.doe.gov/portal/Aaoar1/Aaoar1.download.html'
  remDr$navigate(initial.link)
  remDr$screenshot(TRUE)
  
  # Login protocol
  log_message("Logging into JGI Mycocosm", source_name = source_name)
  
  # Click the login button
  login.button <- '//*[@id="login"]'
  webElem <- remDr$findElement('xpath', login.button)
  webElem$clickElement()
  Sys.sleep(2)
  remDr$screenshot(TRUE)
  
  # Fill the email box
  email.box <- '//*[@id="login"]'
  webElem <- remDr$findElement('xpath', email.box)
  webElem$sendKeysToElement(list(mycocosm_email))
  remDr$screenshot(TRUE)
  
  # Fill the password box
  password.box <- '//*[@id="password"]'
  webElem <- remDr$findElement('xpath', password.box)
  webElem$sendKeysToElement(list(mycocosm_password))
  remDr$screenshot(TRUE)
  
  # Press the sign in button
  sign.in.button <- '//*[@id="home"]/p/input'
  webElem <- remDr$findElement('xpath', sign.in.button)
  webElem$clickElement()
  Sys.sleep(2)
  remDr$screenshot(TRUE)
  
  log_message("Login successful", source_name = source_name)
  
  # Download each genome
  for (k in 1:run.length) {
    run.id <- run.list[k, 'ID']
    run.name <- run.list[k, 'Name']
    go.link <- run.list[k, 'Go.Download.Link']
    
    log_message(paste("Downloading", k, "of", run.length, ":", run.id, "-", run.name), source_name = source_name)
    
    # Create url to Taxa XML directory using the Taxa ID
    xml.link <- str_c(
      'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=',
      run.id,
      '&organizedByFileType=false'
    )
    
    # Navigate to the Taxa XML directory
    Sys.sleep(2)
    remDr$navigate(xml.link)
    Sys.sleep(2)
    
    # Actually download the file (will go to the configured download location)
    Sys.sleep(2)
    remDr$navigate(go.link)
    Sys.sleep(3)
    
    log_message(paste("Downloaded:", run.id), source_name = source_name)
  }
  
  # Close the driver (on.exit will also handle this if function errors)
  tryCatch({
    remDr$close()
  }, error = function(e) {
    # Browser may already be closed, ignore
  })
  
  return()
}

# Download genomes with proper cleanup
log_message(paste("Starting download of", nrow(download_df_subset), "genomes"), source_name = source_name)

# Ensure cleanup happens even if script errors
on.exit({
  tryCatch({
    if (exists("remDr") && !is.null(remDr)) {
      remDr$close()
    }
  }, error = function(e) {
    # Browser may already be closed, ignore
  })
  tryCatch({
    if (exists("rd") && !is.null(rd)) {
      rd$server$stop()
    }
  }, error = function(e) {
    # Server may already be stopped, ignore
  })
}, add = TRUE)

tryCatch({
  download_from_links(download_df_subset, output_dir)
  log_message("Download complete", source_name = source_name)
}, error = function(e) {
  log_message(paste("Error during download:", e$message), source_name = source_name)
  stop(e)
}, finally = {
  # Clean up RSelenium
  tryCatch({
    if (exists("remDr") && !is.null(remDr)) {
      remDr$close()
    }
  }, error = function(e) {
    # Browser may already be closed, ignore
  })
  tryCatch({
    if (exists("rd") && !is.null(rd)) {
      rd$server$stop()
    }
  }, error = function(e) {
    # Server may already be stopped, ignore
  })
})

log_message("Mycocosm download script complete", source_name = source_name)
