# MethMotif Protocol 2024

# Table of Contents

- Preparation
    - [Support Protocol 1](#support-protocol-1)
    - [Support Protocol 2](#support-protocol-2)
    - [Support Protocol 3](#support-protocol-3)
- Analysis
    - [Protocol 1](#protocol-1)
    - [Protocol 2](#protocol-2)
    - [Protocol 3](#protocol-3)
    - [Protocol 4](#protocol-4)

---

## Support Protocol 1

### Software installation

- Install R

    1. To download R, navigate to the Comprehensive R Archive Network (CRAN) [website](https://cran.r-project.org/) and click Download R for Linux/macOS/Windows. Follow the instructions to download the latest base version of R (***4.4.1 as of July 2024).

    2. [OPTIONAL] To install **RStudio**, an integrated development environment (IDE) for R, go to the RStudio [website](https://posit.co/products/open-source/rstudio/) and click on **Download RStudio Desktop**. Then choose the appropriate installer for your operating system and run it. 

- Install the required packages:

    1. Activate R / open RStudio.

    2. Install prerequisites:

        ```r
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")

        BiocManager::install("GenomicRanges")
        BiocManager::install("rGREAT")

        if (!require("remotes", quietly = TRUE))
            install.packages("remotes")
        ```

    3. Install TFregulomeR and ForkedTF:

        ```r
        BiocManager::install("benoukraflab/ TFregulomeR")
        BiocManager::install("benoukraflab/forkedTF")
        ```

- Create directory structure:

    To follow this protocol, create this directory at `~/Projects/tfregulomer_tutorial`. Depending on the local machine's OS, R interprets the `~` key differently.

    - For Windows the full path will be `C:/Users/USERNAME/Documents Projects/tfregulomer_tutorial`.

    - For MacOS the full path will be `/Users/USERNAME/Projects/tfregulomer_tutorial`.

    - For Linux the full path will be `/home/USERNAME/Projects/tfregulomer_tutorial`.

- Download the local database:

    1. Navigate to the [MethMotif.org](https://methmotif.org/API_TFregulomeR/downloads/) and download the latest TFregulomeR database.

        ```bash
        wget https://methmotif.org/API_TFregulomeR/downloads/TFregulomeR_database_2.3.zip
        ```

    2. Unzip the package to a location of your choosing

        ```bash
        mv TFregulomeR_database_2.3.zip ~/Projects/tfregulomer_tutorial/

        cd ~/Projects/tfregulomer_tutorial/

        unzip -q TFregulomeR_database_2.3.zip

        rm TFregulomeR_database_2.3.zip
        ```

    3. Record the file path for the `tfregulomer.sqlite` file, currently it is under this directory `~/Projects/tfregulomer_tutorial/TFregulomeR_database_2.3/`. 

## Support Protocol 2

### Docker installation

- Prepare Docker

    1. To download Docker, navigate to the Get Docker page [here](https://docs.docker.com/get-docker/) and select Docker Desktop for macOS/Windows. Follow the instructions to download the latest Docker Desktop (***4.32.0 as of July 2024).

    2. For Linux users, downloading only the Docker Engine is advised, as Docker Desktop still requires virtualization. Navigate to the docker [docs](https://docs.docker.com/engine/install/) and select the relevant distribution.

- Download TFregulomeR image

    1. Pull the `benoukraflab/tfregulomer` docker image from the Docker Hub repository, using either the Docker Desktop search bar or the command line.

        ```bash
        docker pull benoukraflab/tfregulomer
        ```

- Run TFregulomeR container

    1. Run the container through either Docker desktop.

        ```docker
        docker run -d --mount type=bind,src=/home/USERNAME/Projects,target=/home/rstudio/Projects -e PASSWORD=1234 -p 8787:8787 benoukraflab/tfregulomer
        ```
    
    2. Access the container by opening the browser and navigating to [http://localhost:8787/](http://localhost:8787/).

    3. Login into the RStudio server with the username `rstudio` and the password from the Environment variable.

## Support Protocol 3

### Verify installation

- Setup R environment:

    1. Open the R console from a terminal/command prompt or within RStudio.

    2. Load the TFregulomeR library using the built-in library function.

        ```r
        library(TFregulomeR)
        ```

- Test installation:

    3. Run the dataBrowser function without arguments to check if the package was installed correctly.

        ```r
        all_records <- dataBrowser()
        ```

    4. This should generate a table of the available datasets and their corresponding metadata.

        ```
        2333 record(s) found: ...
        ... covering 676 TF(s)
        ... from 2 species:
        ... ...mouse, human
        ... from 35 organ(s):
        ... ... brain, stem_cell, blood_and_lymph, connective_tissue, liver, colorectum, muscle, bone, stomach, prostate, pancreas, skin, eye, breast, intestine, kidney, lung, esophagus, heart, testis, uterus, spleen, limb, body, cervix, placenta, undefined, adrenal_gland, neck_and_mouth, head, ovary, pleura, thymus, fallopian, vagina
        ... in 3 sample type(s):
        ... ... primary_cells, cell_line, tissue
        ... in  721  different cell(s) or tissue(s)
        ... in 8 type(s) of disease state(s):
        ... ... normal, tumor, Simpson_Golabi_Behmel_syndrome, progeria, metaplasia, unknown, immortalized, premetastatic
        ... from the source(s): GTRD, MethMotif
        ```

    5. If all the steps were followed correctly and TFregulomeR was installed without issue, run the head function to view the first six lines of the dataset:

        ```r
        head(all_records)
        ```

## Protocol 1

### Explore transcription factor dimerization partner

- Setup Rscript

    Here we will load the TFregulomeR library and set up the working directory. By setting the working directory, the results will be stored in a structured and logical manner. The benefit of keeping the local database under the same parent directory is evident when creating a variable to store the full path. With the working directory set, we can create a variable that points to the local database.


    1.	Using a local version of R:

        ```r
        library(TFregulomeR)
        main_dir <- "~/Projects/tfregulome_tutorial"
        dir.create(main_dir)
        db_path <- file.path(main_dir, "TFregulomeR_database_2.3/tfregulome.sqlite")
        setwd(main_dir)
        ```

    1.	Using the Docker version of R:

        ```r
        library(TFregulomeR)
        main_dir <- "~/Projects/tfregulome_tutorial"
        dir.create(main_dir)
        db_path <- Sys.getenv("LOCAL_DATABASE")
        setwd(main_dir)
        ```

- Explore the MethMotif compendium with the TFregulomeR dataBrowser

    The dataBrowser function provides an effective way to query the compendium. The output of the function is a data frame containing relevant metadata. Using `unique` we can print the names of the tissues containing TFBS for CEBPB. The number of available TFs is calculated by filtering the data frame to only include results from the specified cell line and then displaying the number to the console. K562 currently contains the highest number of TF ChIP-seq datasets.


    2.	Examine the extent of CEBPB records in the compendium:

        ```r
        all_records <- dataBrowser(local_db_path = db_path)
        CEBPB_records <- dataBrowser(tf = "CEBPB", species = "human")
        unique(CEBPB_records$cell_tissue_name)
        ```

    3.	Find the tissue with the highest number of available TFs:

        ```r
        for (tissue in CEBPB_records$cell_tissue_name) {
            tf_num <- nrow(all_records[all_records$cell_tissue_name == tissue, ])
            print(paste0(tissue, " has ", tf_num, " TF(s)"))
        }
        ```

- Generate a cofactor list

    To quickly determine the transcription factors that share peak regions, we run a pairwise overlap analysis using `intersectPeakMatrix`. This function computes the union of peaks shared by two transcription factors X and Y along with the overlapping percentage. The percentage of overlapping peaks is computed against the total peaks of X or Y. This percentage is first determined by which transcription factor dataset is being accessed. If there are 500 shared peaks, and TF X has 10,000 peaks, the overlap amounts to 5% of TF X's peaks. However, if TF Y has 5000 peaks in total, then the overlap amount to 10% of Y's peak. The same is true for MethMotif logos created from the X∩Y peaks; the binding motif is derived from the TF X overlapping peaks. Similarly, the logo generated from the Y∩X overlapping peaks is derived from TF Y overlapping peaks.

    When running `intersectPeakMatrix` on TF lists containing over ten factors, using the `local_db_path` option should decrease the execution time by over 10x. Due to the long run time required to compute intersect matrices, it is advised to save the variable to the local machine, using R's built-in function `saveRDS`. Hence, when running this analysis again into the same session, the stored matrix can be read without rerunning Step 6.


    4.	Create a subdirectory for storing the output of the K562 CEBPB analysis:

        ```r
        sub_dir <- "k562_cebpb"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        ```

    5.	Explore the complete list of transcription factors and the metadata within the TFregulomeR compendium belonging to the K562 cell line:

        ```r
        k562_TFBS <- dataBrowser(cell_tissue_name = "K562")
        unique(k562_TFBS$TF)
        ```

    6.	Investigate the overlapping peak regions between CEBPB and all other TFs in the K562 cell line:

        ```r
        intersectMatrix_k562_CEBPB <- intersectPeakMatrix(
            peak_id_x = c("MM1_HSA_K562_CEBPB"),
            motif_only_for_id_x = TRUE,
            peak_id_y = k562_TFBS$ID,
            motif_only_for_id_y = TRUE,
            local_db_path = db_path
        )
        ```

    7.	Save the computed intersect matrix:

        ```r
        saveRDS(intersectMatrix_k562_CEBPB, file = "intersectMa-trix_k562_cebpb.rds")
        ```

- Extract cofactor information

    Given the two-dimensional structure of the output, a function is required to process the results in a given direction. In this example, CEBPB peaks overlaps with peaks of all other TFs that are stored in the X direction, resulting in a slice with dimensions 1 x 131. While in the Y direction, a slice with dimensions 131 x 1 points to all the TFs analyzed in K562 with overlaps with CEBPB. The original output being in alphabetical order, it is necessary to rearrange the data according to the co-binding percentage. This can be achieved by reformatting the intersection matrix into a data frame, setting it in descending order. Printing the head of this ordered data frame returns the names and co-binding percentages. The result function further provides a new MethMotif logo for each TF available in the selected direction. In this example, ten PDF files were generated and stored in the current directory.

    If RStudio was installed, the plots are available in the bottom right panel under "Plots." Using the arrow tabs to navigate through the plots, it is easy to see the subtle differences in each logo.


    8.	Read the matrix in the X direction to understand CEBPB partners:
        
        ```r
        k562_CEBPB_result <- intersectPeakMatrixResult(
            intersectPeakMatrix = intersectMatrix_k562_CEBPB,
            return_intersection_matrix = TRUE,
            angle_of_matrix = "x"
        )  
        ```

    9.	Determine the top 10 most common partners of CEBPB in K562:

        ```r
        k562_CEBPB_result_t <- as.data.frame(t(k562_CEBPB_result$intersection_matrix))
        attach(k562_CEBPB_result_t)
        k562_CEBPB_result_order <- as.data.frame(
            k562_CEBPB_result_t[order(-MM1_HSA_K562_CEBPB),, drop = FALSE]
        )
        detach(k562_CEBPB_result_t)
        head(k562_CEBPB_result_order, n = 10)
        ```

    10.	Determine the top 10 most common partners of CEBPB in K562:

        ```r
        top_cobinding <- rownames(head(k562_CEBPB_result_order, n = 10))
        intersectPeakMatrixResult(
            intersectPeakMatrix = intersectMatrix_k562_CEBPB,
            save_MethMotif_logo = TRUE,
            angle_of_logo = "x",
            saving_MethMotif_logo_y_id = top_cobinding
        )
        ```

## Protocol 2

### Alternative cofactor visualization

- Report based visualization

    `miniCofactorReport` uses the ReMapEnrich (Chèneby *et al.*, 2020) functionality to calculate an adjusted p-value for each subset of peaks. Following the statistical test, it updates the cofactor ranking based on the new score.

    1.	To better visualize bZIP co-binding load the forkedTF library:

        ```r
        library(forkedTF)
        ```

    2.	The forkedTF library contains a report function that summarizes cofactor information:

        ```r
        miniCofactorReport(
            TF = "CEBPB",
            cell = "K562",
            filterBy = "q.significance",
            Methylation = TRUE,
            includeMotifOnly = TRUE,
            pdfName = "CEBPB_K562_miniCofactorReport.pdf",
            local_db_path = db_path
        )
        ```

- FPWM visualization

    We can use the R function `lapply` to build a list of TFs by extracting the last part of the TFregulomeR IDs (MM1_HSA_K562_CEBPB), from the variable `top_cobinding` generated in Step 10 above. The `createFPWM` function builds a position weight matrix, allowing for a portion of the matrix to be calculated in independent blocks based on a subset of the peaks. This results in a consistent TF matrix half forking to several partner-specific half matrices. This representation thereby deconvolutes the degenerate section of bZIP protein motifs. Finally, we construct a graphical representation of the forked position weight matrix. The representation excels at segregating TF partner motifs compared to the primary TF motif.

    3.	Build a list of transcription factor names to fork against CEBPB:

        ```r
        tfs <- lapply(top_cobinding[-1], function(x) {
            return(strsplit(x, "_")[[1]][4])
        })
        ```

    4.	From the forkedTF library, execute the createFPWM function:

        ```r
        fpwm <- createFPWM(
            mainTF = "CEBPB",
            partners = tfs,
            cell = "K562",
            forkPosition = 5,
            flipMatrix = FALSE,
            local_db_path = db_path
        )
        ```

    5.	Plot the PWM generate and provide a file name for the PDF:

        ```r
        plotFPWM(fpwm, pdfName = "MM1_HSA_K562_CEBPB_fpwm_plot.pdf")
        write.FPWM(
            FPWM = fpwm,
            format = "transfac",
            fileName = "MM1_HSA_K562_CEBPB_fpwm.tf"
        )
        ```

## Protocol 3

### Characterize bZIP Partners/Co-Factors

- Prepare cofactor specific analysis

    We will create a directory under "k562_cebpb" called "overlapped_with_ATF4". This subdirectory will hold the output of the characterization of TFBS, where both CEBPB and ATF4 are present.


    1. Create a subdirectory to store the MethMotif logos from CEBPB binding sites that overlap with ATF4 binding:

        ```r
        sub_dir <- "k562_cebpb/overlapped_with_ATF4"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        ```

- Cofactor ATF4

    The function `commonPeaks` computes the locations shared between a reference TF and a list of target TFs. Only the peak regions shared between the listed TFs will be output. To reduce complexity when handling the output from `commonPeaks`, use the `commonPeaksResult` function. This allows for granular control over what information is returned.  `commonPeaksResult` only reformats the output of `commonPeaks` and therefore doesn't generate any new data. For example, returning the methylation profile from `commonPeaksResult` will only succeed if the narrow methylation profile was calculated in the preceding execution of `commonPeaks`.

    2. Select all the DNA binding sites where CEBPB and ATF4 peaks overlapped:

        ```r
        k562_CEBPB_com <- commonPeaks(
            target_peak_id = "MM1_HSA_K562_CEBPB",
            motif_only_for_target_peak = TRUE,
            compared_peak_id = c("MM1_HSA_K562_ATF4"),
            motif_only_for_compared_peak = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        ```

    3. Extract information from subset of overlapping peaks:

        ```r
        k562_CEBPB_com_res <- commonPeakResult(
            commonPeaks = k562_CEBPB_com,
            return_common_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE,
            return_summary = TRUE
        )
        ```

    4. Store the results into new variables:

        ```r
        peak_summary <- k562_CEBPB_com_res$peak_summary
        common_peak_list <- k562_CEBPB_com_res$common_peak_list
        methylation_profile <- k562_CEBPB_com_res$methylation_profile
        ```

    5.	Visualize the percentage of CEBPB peaks where ATF4 dimerization occurs:

        ```r
        peak_summary
        ```

    6.	Preview the list of common peak locations: 

        ```r
        K562_CEBPB_com_peaks <- common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
        head(K562_CEBPB_com_peaks)
        ```

    7.	View the narrow range beta scores in increments of 10%: 

        ```r
        methylation_profile$MM1_HSA_K562_CEBPB_common_peaks
        ```

- Annotate CEBPB∩ATF4

    In order to facilitate downstream analysis of the results generated by the `commonPeaks` function, we can use the `exportMMPFM` to save the common peak regions as position frequency matrices. TFregulomeR offers basic analysis tools, such as the rGREAT package, to compute GO term enrichment.

    8.	Save a copy of the Frequency Position and Beta Score matrices:

        ```r
        exportMMPFM(
            fun_output = k562_CEBPB_com,
            fun = "commonPeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        ```
   
    9.	Annotate the common peaks with genomic location information, i.e., promoter-TSS, intron, exon, etc: 

        ```r
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        k562_CEBPB_com_peaks_loc <- genomeAnnotate(
            peaks = K562_CEBPB_com_peaks,
            return_annotation = TRUE,
            return_html_report = TRUE,
            local_db_path = db_path
        )
        ```
    
    10.	Determine the enriched Gene Ontology terms for the common peaks: 

        ```r
        library(rGREAT)
        library(plotly)
        library(crosstalk)
        library(DT)
        K562_CEBPB_com_peaks_func <- greatAnnotate(
            peaks = K562_CEBPB_com_peaks,
            return_annotation = TRUE,
            return_html_report = TRUE
        )
        ```

- Specific analysis for CEBPB/CEBPD heterodimer

    Reiterate previous steps 1 to 10 for the CEBPB/CEBPD heterodimer.

    11.	Read the matrix in the X direction to understand CEBPB partners:

        ```r
        sub_dir <- "k562_cebpb/overlapped_with_CEBPD"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        k562_CEBP_B_D_com <- commonPeaks(
            target_peak_id = "MM1_HSA_K562_CEBPB",
            motif_only_for_target_peak = TRUE,
            compared_peak_id = c("MM1_HSA_K562_CEBPD"),
            motif_only_for_compared_peak = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        k562_CEBP_B_D_com_res <- commonPeakResult(
            commonPeaks = k562_CEBP_B_D_com,
            return_common_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE,
            return_summary = TRUE
        )
        peak_summary <- k562_CEBP_B_D_com_res$peak_summary
        common_peak_list <- k562_CEBP_B_D_com_res$common_peak_list
        methylation_profile <- k562_CEBP_B_D_com_res$methylation_profile
        peak_summary
        K562_CEBP_B_D_com_peaks <- common_peak_list$MM1_HSA_K562_CEBPB_common_peaks
        head(K562_CEBP_B_D_com_peaks)
        methylation_profile$MM1_HSA_K562_CEBPB_common_peaks
        exportMMPFM(
            fun_output = k562_CEBP_B_D_com,
            fun = "commonPeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        k562_CEBP_B_D_com_peaks_loc <- genomeAnnotate(
            peaks = K562_CEBP_B_D_com_peaks,
            return_annotation = TRUE,
            return_html_report = TRUE,
            local_db_path = db_path
        )
        K562_CEBP_B_D_com_peaks_func <- greatAnnotate(
            peaks = K562_CEBP_B_D_com_peaks,
            return_annotation = TRUE,
            return_html_report = TRUE
        )
        ```

- Global analysis of genomic locations of CEBPB peak in K562 cells

    Using the `loadPeaks` function, we can access the original K562 CEBPB peaks in TFregulomeR. Again, ensure only to **generate the HTML report in a specific subdirectory**, to avoid overwriting previous reports.

    12.	Before analyzing the location and function information, create a baseline for all CEBPB peaks:

        ```r
        sub_dir <- "k562_cebpb"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        k562_CEBPB_peaks <- loadPeaks(
            id = "MM1_HSA_K562_CEBPB",
            includeMotifOnly = TRUE,
            local_db_path = db_path
        )
        k562_CEBPB_peaks_loc <- genomeAnnotate(
            peaks = k562_CEBPB_peaks,
            return_annotation = TRUE,
            return_html_report = TRUE,
            local_db_path = db_path
        )
        K562_CEBPB_peaks_func <- greatAnnotate(
            peaks = k562_CEBPB_peaks,
            return_annotation = TRUE,
            return_html_report = TRUE
        )
        ```

- Analyze the differences between the CEBPB dimers

    From the three donut plots embedded in the HTML files, CEBPB dimerized with ATF4 generally follows the global trend of primarily intronic and, to a lesser extent, the intergenic regions. The result is the plot for CEBPB dimerized to CEBPD, where the primary location is intronic regions, and the secondary site is promoter-TSS regions. Interestingly, CEBPB-ATF4 had about 10% of its TFBS methylated, while CEBPB-CEBPD had 80% of its TFBS methylated. The relevance of binding methylated DNA in the promoter-TSS region is "specific pioneer TFs in reshaping the chromatin landscape in cancer patients to rewire gene regulatory networks through local DNA demethylation of cis-regulatory regions" (Lemma *et al.*, 2022).    This information will be built on in the next step. From the enriched terms and their corresponding adjusted p-value, it is possible to determine what genes each specific CEBPB dimer is targeting. Both the baseline and CEBPB-ATF4 peaks share most of their top terms. However, the CEBPB-CEBPD TFBS has "nucleic acid metabolic process" and "nucleobase-containing compound metabolic process" as the leading terms.

## Protocol 4

### Context-independent and dependent analysis

- Prepare MAFF analysis

    In using `dataBrowser` we extract a single record from every tissue that contains MAFF. Currently, this returns three records. Then using the ID metadata from the `dataBrowser` output as the input for `searchMotif` returns a custom R class. This custom R class contains the relevant information to plot a MethMotif logo. With this, a baseline plot for each cell line is generated.

    1. Create a subdirectory to store the context-independent and dependent analysis for MAFF binding sites:

        ```r
        sub_dir <- "maff_example"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        ```

    2. Query the TFrefulomeR compendium for occurrences of the transcription factor MAFF in human cell lines and tissues:

        ```r
        MAFF_record <- dataBrowser(
            tf = "MAFF",
            species = "human",
            local_db_path = db_path
        )
        ```

    3. Generate a unique MethMotif logo for each cell line that contains MAFF: 
        
        ```r
        for (tf_id in MAFF_record$ID) {
            tf_motif <- searchMotif(id = tf_id)
            plotLogo(MM_object = tf_motif)
        }
        ```

- Browse the plots in RStudio or open them from the file system:

    There are two contexts to focus on in the MAFF MethMotif logos. First is an enriched TGA motif before the primary MAFF motif found only in the K562 cell line.  The second context is the methylation level within the motifs' first 6-7 bases. The K562 motif methylation levels are lower than those computed for the Hela-s3 and hepG2 cell lines. This indicates that the partner transcription factor responsible for the TGA motif selectively binds unmethylated DNA.

- Explore context independent MAFF

    The peak regions shown to be shared among all cell lines can be considered the context-independent targets of MAFF. The motifs generated from the shared peaks may vary because the "X" and "Y" directions in the TF intersection matrix are not identical. The peak set is shared, but the MethMotif logo is built from the corresponding ChIP-seq experiment for each cell line. The peak summary provides the percentage of overlap between the common peaks and the MAFF peaks for each cell line. Within the `peak_list` and `methylation_profile`, each cell line has a data frame. This data frame is extracted using the `$` operator followed by the TFregulomeR ID.

    4. Create a subdirectory to store the context-independent analysis for MAFF binding sites:

        ```r
        sub_dir <- "maff_example/common_peaks"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        ```

    5. Discover the shared peak regions between all MAFF TFBS:

        ```r
        all_MAFF_com <- commonPeaks(
            target_peak_id = MAFF_record$ID,
            motif_only_for_target_peak = TRUE,
            compared_peak_id = MAFF_record$ID,
            motif_only_for_compared_peak = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        ```

    6. Create a MethMotif logo for each cell line:

        ```r
        all_MAFF_com_res <- commonPeakResult(
            commonPeaks = all_MAFF_com,
            return_common_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE,
            return_summary = TRUE
        )
        ```

    7. Categorize the results into new variables:

        ```r
        all_MAFF_com_peak_summary <- all_MAFF_com_res$peak_summary
        all_MAFF_com_peak_list <- all_MAFF_com_res$common_peak_list
        all_MAFF_com_methylation_profile <- all_MAFF_com_res$methylation_profile
        ```

    8. Review the common peak results:

        ```r
        all_MAFF_com_peak_summary
        HeLa_MAFF_common_peaks <- all_MAFF_com_peak_list$`MM1_HSA_HeLa-S3_MAFF_common_peaks`
        HepG2_MAFF_common_peaks <- all_MAFF_com_peak_list$MM1_HSA_HepG2_MAFF_common_peaks
        K562_MAFF_common_peaks <- all_MAFF_com_peak_list$MM1_HSA_K562_MAFF_common_peaks
        all_MAFF_com_methylation_profile$`MM1_HSA_HeLa-S3_MAFF_common_peaks`
        all_MAFF_com_methylation_profile$MM1_HSA_HepG2_MAFF_common_peaks
        all_MAFF_com_methylation_profile$MM1_HSA_K562_MAFF_common_peaks
        ```

    9. Save a copy of the Frequency Position and Beta Score matrices:

        ```r
        exportMMPFM(
            fun_output = all_MAFF_com,
            fun = "commonPeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        ```

- Context independent MAFF peaks in K562

    We can now compare the subset peaks that are shared across all cell lines within the K562 cell line. The cofactor report provides the most relevant information from the `intersectPeakMatrix` analysis regarding the overlapped regions: percentage of overlap, motif, methylation level, and read enrichment for each cofactor in descending order. The motifs built for HeLa-s3 and hepG2 using the common peak regions look nearly identical. Compare this to the motif created in K562 using the same common peak regions. There is less hypermethylation and more hypomethylation. Besides methylation, the motifs are all but identical. This is expected for the context-independent binding sites of a transcription factor.

    10. Investigate the cofactors found in common TFBS with the K562 cell line:

        ```r
        K562_MAFF_com_peaks_intersect <- intersectPeakMatrix(
            user_peak_list_x = list(K562_MAFF_common_peaks),
            user_peak_x_id = "MM1_HSA_K562_MAFF",
            peak_id_y = k562_TFBS$ID,
            motif_only_for_id_y = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        saveRDS(
            K562_MAFF_com_peaks_intersect,
            file = "K562_maff_com_peaks_intersect.rds"
        )
        ```

    11. Build a report summarizing all the cofactors from the intersectPeakMatrix:

        ```r
        cofactorReport(
            intersectPeakMatrix = K562_MAFF_com_peaks_intersect,
            cobinding_threshold = 0.1
        )
        ```

- Context dependent MAFF peaks in K562

    The function `exclusivePeaks` helps isolate biological contexts between a target TF and a list of TFs to compare against. Only the peak regions unique to the target TF will be output. Leveraging this function to remove the peaks shared among MAFF in all cell lines leaves only the peaks specific to K562. This function forms the backbone of context-dependent analysis. The variable `K562_MAFF_exclu_peaks` will contain the K562-specific peaks for MAFF in a BED-like format. To facilitate downstream analysis of the results generated by the `exclusivePeaks` function, use the `exportMMPFM` function.

    12. Create a subdirectory to store the context-dependent analysis for MAFF binding sites in K562:

        ```r
        sub_dir <- "maff_example/exclusive_peaks/k562"
        dir.create(file.path(main_dir, sub_dir), recursive = TRUE)
        setwd(file.path(main_dir, sub_dir))
        ```

    13. Consider the peak regions unique to MAFF TFBS within K562:

        ```r
        K562_MAFF_exclu <- exclusivePeaks(
            target_peak_id = "MM1_HSA_K562_MAFF",
            motif_only_for_target_peak = TRUE,
            excluded_peak_id = c("MM1_HSA_HeLa-S3_MAFF", "MM1_HSA_HepG2_MAFF"),
            motif_only_for_excluded_peak = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        ```

    14. Separate the DNA binding sites exclusive to MAFF in K562 and plot a MethMotif logo for these regions:

        ```r
        K562_MAFF_exclu_res <- exclusivePeakResult(
            exclusivePeaks = K562_MAFF_exclu,
            return_exclusive_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE
        )
        K562_MAFF_exclu_peaks <- K562_MAFF_exclu_res$exclusive_peak_list$MM1_HSA_K562_MAFF_exclusive_peaks
        ```

    15. Save a copy of the Frequency Position and Beta Score matrices:

        ```r
        exportMMPFM(
            fun_output = K562_MAFF_exclu,
            fun = "exclusivePeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        ```

    16. Investigate the cofactors that are exclusive to MAFF TFBS within the K562 cell line:

        ```r
        K562_MAFF_exclu_peaks_intersect <- intersectPeakMatrix(
            user_peak_list_x = list(K562_MAFF_exclu_peaks),
            user_peak_x_id = "MM1_HSA_K562_MAFF",
            peak_id_y = k562_TFBS$ID,
            motif_only_for_id_y = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        saveRDS(
            K562_MAFF_exclu_peaks_intersect,
            file = "K562_maff_exclu_peaks_intersect.rds"
        )
        ```

    17. Display the cofactor information for K562-specific MAFF binding sites:

        ```r
        cofactorReport(
            intersectPeakMatrix = K562_MAFF_exclu_peaks_intersect,
            cobinding_threshold = 0.1
        )
        ```

- Repeat exclusive peak analysis

    To facilitate the comparison of context independent and dependent peaks we can repeat the above steps for the other two cell lines.

    18. Analyze the exclusive peak regions of MAFF in HeLa-s3 and characterize cell-specific cofactors:

        ```r
        sub_dir <- "maff_example/exclusive_peaks/hela-s3"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        HeLa_MAFF_exclu <- exclusivePeaks(
            target_peak_id = "MM1_HSA_HeLa-S3_MAFF",
            motif_only_for_target_peak = TRUE,
            excluded_peak_id = c("MM1_HSA_K562_MAFF", "MM1_HSA_HepG2_MAFF"),
            motif_only_for_excluded_peak = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        HeLa_MAFF_exclu_res <- exclusivePeakResult(
            exclusivePeaks = HeLa_MAFF_exclu,
            return_exclusive_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE
        )
        HeLa_MAFF_exclu_peaks <- HeLa_MAFF_exclu_res$exclusive_peak_list$`MM1_HSA_HeLa-S3_MAFF_exclusive_peaks`
        exportMMPFM(
            fun_output = HeLa_MAFF_exclu,
            fun = "exclusivePeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        HeLa_TFBS <- dataBrowser(cell_tissue_name = "HeLa-S3")
        HeLa_MAFF_exclu_peaks_intersect <- intersectPeakMatrix(
            user_peak_list_x = list(HeLa_MAFF_exclu_peaks),
            user_peak_x_id = "MM1_HSA_HeLa-S3_MAFF",
            peak_id_y = HeLa_TFBS$ID,
            motif_only_for_id_y = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        saveRDS(
            HeLa_MAFF_exclu_peaks_intersect,
            file = "hela_maff_exclu_peaks_intersect.rds"
        )
        cofactorReport(
            intersectPeakMatrix = HeLa_MAFF_exclu_peaks_intersect,
            cobinding_threshold = 0.1
        )
        ```

    19. Analyze the exclusive peak regions of MAFF in HepG2 and characterize cell-specific cofactors:

        ```r
        sub_dir <- "maff_example/exclusive_peaks/hepg2"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        HepG2_MAFF_exclu <- exclusivePeaks(
            target_peak_id = "MM1_HSA_HepG2_MAFF",
            motif_only_for_target_peak = TRUE,
            excluded_peak_id = c("MM1_HSA_K562_MAFF", "MM1_HSA_HeLa-S3_MAFF"),
            motif_only_for_excluded_peak = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        HepG2_MAFF_exclu_res <- exclusivePeakResult(
            exclusivePeaks = HepG2_MAFF_exclu,
            return_exclusive_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE
        )
        HepG2_MAFF_exclu_peaks <- HepG2_MAFF_exclu_res$exclusive_peak_list$MM1_HSA_HepG2_MAFF_exclusive_peaks
        exportMMPFM(
            fun_output = HepG2_MAFF_exclu,
            fun = "exclusivePeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        HepG2_TFBS <- dataBrowser(cell_tissue_name = "HepG2")
        HepG2_MAFF_exclu_peaks_intersect <- intersectPeakMatrix(
            user_peak_list_x = list(HepG2_MAFF_exclu_peaks),
            user_peak_x_id = "MM1_HSA_HepG2_MAFF",
            peak_id_y = HepG2_TFBS$ID,
            motif_only_for_id_y = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        saveRDS(
            HepG2_MAFF_exclu_peaks_intersect,
            file = "hepg2_maff_exclu_peaks_intersect.rds"
        )
        cofactorReport(
            intersectPeakMatrix = HepG2_MAFF_exclu_peaks_intersect,
            cobinding_threshold = 0.1
        )
        ```

- Review the cofactor reports generated

    The HepG2 cofactor report contains only one cofactor, as MAFF∩MAFF is the same as the global motif. This is due to the high number of peaks preventing other cofactors from reaching the 10% overlap threshold. Within the HeLa-s3 cofactor report, the cofactors JUN, JUND, and FOS show slight enrichment for the TGA motif. However, NFE2L2 displays the most robust enrichment of the TGA motif. Notice the change in methylation and its association with the TGA motif. Finally, the K562 cofactor report displays a strong TGA enrichment in all cofactor motifs. Since NFE transcription factors contain the motif in two cell lines, the next step is to remove peaks shared between MAFF and NFEs.

- MAFF in K562 without NFE factors

    This analysis aims to further subset the K562-specific MAFF binding sites by removing any NFE transcription factors. If NFE cobinding was responsible for the TGA motif, the MethMotif logo saved should not have TGA enrichment at sites 3-5. This is a subset of a subset of the total peak region for MAFF within K562. By removing two of the main cofactors from the exclusive peak regions, it is expected that the order and the motifs of the remaining cofactors will have changed.

    20. Create a subdirectory to store the context-dependent analysis for MAFF binding sites in K562 without NFE factors:

        ```r
        sub_dir <- "maff_example/exclusive_peaks/k562/no_nfe"
        dir.create(file.path(main_dir, sub_dir))
        setwd(file.path(main_dir, sub_dir))
        ```

    21. Remove the NFE TFs peak regions from the K562-specific MAFF TFBS:

        ```r
        K562_MAFF_exclu_no_NFE <- exclusivePeaks(
            user_target_peak_list = list(K562_MAFF_exclu_peaks),
            user_target_peak_id = "MM1_HSA_K562_MAFF",
            excluded_peak_id = c("MM1_HSA_K562_NFE2", "MM1_HSA_K562_NFE2L2"),
            motif_only_for_excluded_peak = TRUE,
            local_db_path = db_path
        )
        K562_MAFF_exclu_no_NFE_res <- exclusivePeakResult(
            exclusivePeaks = K562_MAFF_exclu_no_NFE,
            return_exclusive_peak_sites = TRUE,
            save_MethMotif_logo = TRUE,
            return_methylation_profile = TRUE
        )
        K562_MAFF_exclu_no_NFE_peaks <- K562_MAFF_exclu_no_NFE_res$exclusive_peak_list$MM1_HSA_K562_MAFF_exclusive_peaks
        ```

    22. Save a copy of the Frequency Position and Beta Score matrices:

        ```r
        exportMMPFM(
            fun_output = K562_MAFF_exclu_no_NFE,
            fun = "exclusivePeaks",
            save_motif_PFM = TRUE,
            save_betaScore_matrix = TRUE
        )
        ```

    23. Test the cofactors within the no NFE exclusive peaks:

        ```r
        K562_MAFF_exclu_no_NFE_cofactor <- intersectPeakMatrix(
            user_peak_list_x = list(K562_MAFF_exclu_no_NFE_peaks),
            user_peak_x_id = "MM1_HSA_K562_MAFF",
            peak_id_y = k562_TFBS$ID,
            motif_only_for_id_y = TRUE,
            methylation_profile_in_narrow_region = TRUE,
            local_db_path = db_path
        )
        saveRDS(
            K562_MAFF_exclu_no_NFE_cofactor,
            file = "K562_maff_exclu_no_nfe_cofactor.rds"
        )
        ```

    24. Build the new cofactor report on the final peak set:

        ```r
        cofactorReport(K562_MAFF_exclu_no_NFE_cofactor, cobinding_threshold = 0.1)
        ```

