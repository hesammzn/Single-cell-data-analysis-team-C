#Loads all necessary libraries
library(shiny)
library(bslib)
library(DT)
library(scatterplot3d)

# Loads all data 
tsne_data <- read.csv("data/tsne_3d_data.csv")
number_cells <- read.csv("data/Number_cells.csv", header = FALSE) #removes false headers
colnames(number_cells) <- NULL #removes names for columns 
top20_old <- read.csv("data/Top20_Old.csv")
top_markers_alt <- read.csv("data/top_markers_per_cluster_ALT.csv")
scatterplot_data <- read.csv("data/scatterplot3d_with_colors.csv")
last_table <- read.csv("data/last_table.csv")
colnames(last_table) <- gsub("\\.", " ", colnames(last_table)) #Removes any "." that are found in the headers


# Defines the apps theme and font by using Bootstrap
my_theme <- bs_theme(
  version = 5,
  bg = "#1E1E2F",
  fg = "#FFFFFF",
  primary = "#00AEEF",
  base_font = font_google("Poppins", local = FALSE),
  heading_font = font_google("Montserrat", local = FALSE)
)

###UI Section###
#Creates a navigation bar
ui <- page_navbar(
  id = "navbar",
  title = div(class = "navbar-brand fw-bold", "Team C Website"), #title of the navigation bar
  theme = my_theme, #applies the bootstrap theme
  #Customizes the CSS styles e.g colour of navigation bar or the sizes of cards and panels that were used or colours of the action buttons etc. 
  tags$head(
    tags$style(HTML("
      html, body {
        background-color: #1E1E2F;
      }
      .navbar {
        background-color: #121212 !important;
        border-bottom: 1px solid #4A4A4A;
      }
      .nav-link, .navbar-brand {
        color: #FFFFFF !important;
      }
      .nav-link.active {
        color: #00AEEF !important;
      }
      .dataTables_wrapper,
      table.dataTable {
        background-color: #2C2C3A;
        color: #FFFFFF;
      }
      .dataTables_filter input {
        background-color: #121212;
        color: #FFFFFF;
        border: 1px solid #555;
        border-radius: 5px;
      }
      .dataTables_paginate .paginate_button {
        background-color: #121212;
        color: #FFFFFF !important;
        border: none;
      }
      .dataTables_info {
        color: #AAAAAA;
      }
      .dataTables_filter {
        float: left !important;
        text-align: left !important;
      }

      .big-card {
        width: 300px;
        height: 150px;
        padding: 20px;
        border-radius: 20px;
        display: flex;
        justify-content: center;
        align-items: center;
        text-align: center;
        font-size: 1.5rem;
        margin: 15px;
        transition: background-color 0.3s, transform 0.3s;
      }
      .big-card:hover {
        transform: translateY(-5px);
      }
      .btn-primary {
        background-color: #00AEEF;
        border: none;
      }
      .btn-primary:hover {
        background-color: #008CBA;
      }
      .btn-danger {
        background-color: #FF1919;
        border: none;
      }
      .btn-danger:hover {
        background-color: #CC0000;
      }
      .btn-purple {
        background-color: #B400F1;
        border: none;
      }
      .btn-purple:hover {
        background-color: #8A00C8;
      }
    "))
  ),

###Home page###
#Creates a seperate page for Home and decides its dimensions
  nav_panel("Home",
            tags$div(
              class = "d-flex flex-column align-items-center justify-content-center",
              style = "min-height: 80vh;",
              tags$main(
                class = "px-3 text-center",
#Creates a headline and decides its gap to the next section of the page. Also customises the colour of the font.                
                tags$h1(
                  style = "margin-bottom: 80px;", 
                  "Re-implementation of single-cell RNA-seq pipeline in ",
                  tags$a(
                    href = "https://www.pnas.org/doi/full/10.1073/pnas.1507125112", #a link to the research paper
                    target = "_blank",
                    style = "color: #00AEEF; text-decoration: underline;",
                    "Darmanis et al. (2015)"
                  ),
                  " paper as well as using a newer alternative pipeline published after 2016"
                ),
#Creates a paragraph with more of the same customisations                  
                tags$p(
                  class = "lead",
                  "In this website, we show our findings for analysis on the single-cell RNA-sequencing data from the Darmanis et al. (2015) paper using two pipelines:"
                ),
                
                tags$div(
                  style = "margin-bottom: 60px;",
                  tags$p(
                    style = "margin: 5px; font-size: 1.3rem;",  
                    "1 - Original Analysis: Provided in the supplementary information of the paper."
                  ),
                  tags$p(
                    style = "margin: 5px; font-size: 1.3rem;", 
                    "2 - Alternative Analysis: Using newer tools published after 2016."
                  )
                ),
                
                tags$p(class = "lead", "Use the buttons below to explore the different sections."),
#Creates action buttons which if clicked on will link to the mentioned page.                 
                tags$div(
                  class = "d-flex justify-content-center flex-wrap",
                  actionButton(
                    "go_original", 
                    div(style = "font-weight: bold;", "Original Analysis"), #will link to the original analysis page
                    class = "big-card btn btn-primary"
                  ),
                  actionButton(
                    "go_alternative", 
                    div(style = "font-weight: bold;", "Alternative Analysis"), #will link to alternative analysis page
                    class = "big-card btn btn-danger"
                  ),
                  actionButton(
                    "go_contact", 
                    div(style = "font-weight: bold;", "Contact"), #will link to contacts page
                    class = "big-card btn btn-purple"
                  )
                )
              )  
            )
  ),
  
#Creates another page for the Original analysis  
  nav_panel("Original Analysis",
            card("Re-implementation of pipeline provided in Darmanis et al. (2015):",
                 
                 tags$h4("Methods and Materials"), #header for methods
#Creates a list of sentences for the methods section, that will each be on a separate line. Many of these are linked with a site that references the tool that was used for that part of the analysis.                  
                 tags$ul(
                   tags$li(HTML('Downloading FASTQ files using <a href="https://github.com/ncbi/sra-tools" target="_blank" style="color:#00AEEF;">SRA toolkit</a> (v.3.0.0)')),
                   tags$li(HTML('Trimming sequences with <a href="https://prinseq.sourceforge.net/" target="_blank" style="color:#00AEEF;">Prinseq</a> (v.0.20.4)')),
                   tags$li(HTML('Removing Nextera adapters with <a href="https://github.com/FelixKrueger/TrimGalore" target="_blank" style="color:#00AEEF;">Trim Galore</a> (v.0.6.10)')),
                   tags$li(HTML('Creating genome indices with <a href="https://github.com/alexdobin/STAR/tree/master" target="_blank" style="color:#00AEEF;">STAR</a> (v2.7.11a) for hg19 and GRCh37.p13')),
                   tags$li(HTML('Aligning reads to hg19 using STAR</a>')),
                   tags$li("Merging gene count indices and replacing gene_id with gene_names"),
                   tags$li(HTML('Calculating pairwise distance between cells with <a href="https://www.bioconductor.org/packages/release/bioc/html/scde.html" target="_blank" style="color:#00AEEF;">SCDE</a> (v.2.27.1)')),
                   tags$li(HTML('Dimension reduction with <a href="http://factominer.free.fr/" target="_blank" style="color:#00AEEF;">FactoMineR</a> (v.2.11) and <a href="https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html" target="_blank" style="color:#00AEEF;">t-SNE</a> (v.0.1.3.1)')),
                   tags$li(HTML('Clustering with <a href="https://mclust-org.github.io/mclust/" target="_blank" style="color:#00AEEF;">mclust</a> (v.6.0.0)')),
                   tags$li("Determining top 20 genes expressed in each cluster")
                 ),
                 
                 tags$div(
                   style = "margin-top: 40px; ", #gap to the above section
                   tags$h4("Challenges"), #header for challenges
#Creates list of Challenges from this section of the analysis. 
                   tags$ul(
                     tags$li("No version was mentioned for any of the tools and packages, which can bias the results during re-implementation."),
                     tags$li("The scde package was not running initially, and we had to downgrade some libraries to make scde work."),
                     tags$li("Running scde.error.model was not possible using RStudio, so we had to write SLURM scripts for this part."),
                     tags$li("We did not know how they cleaned the initial gene matrix, including which function and settings were used."),
                     tags$li("No parameters were mentioned for packages during the data analysis part."),
                     tags$li("The hg19 genome has been receiving many updates (such as lnc-RNAs), which can cause differences in results."),
                     tags$li("No details were mentioned about the actual steps of the analysis in the paper; only the packages and libraries used were listed."),
                     tags$li("Some packages use random numbers (such as t-SNE), which can produce different findings each time.")
                   )
                 ),
                 
               
                 tags$div(
                   style = "margin-top: 40px;",
                   tags$h4("Results"), #header 
                   tags$p(
                     "We could not identify clusters since there were many similar genes among the top 20 expressed genes in various clusters. Possible reasons are mentioned above. Other findings of the paper could not be reproduced either due to biased and wrong results from data analysis."
                   )
                 ),
                 
                 
                 # Table 1
                 tags$div(style = "margin-top: 30px; text-align: center;",
                          tags$h5("Table 1. number of cells in each of clusters among different datasets based on original analysis."),
                          dataTableOutput("number_cells_table"),
                 ),
                 
                 # Table 2
                 tags$div(style = "margin-top: 30px; text-align: center;",
                          tags$h5("Table 2. The top 20 expressed genes in each cluster for the gene matrix cleaned with excluding rows with less than 100 counts analyzed with pipeline provided in the paper."),
                          dataTableOutput("top20_old_table"),
                 ),
                 
                 # Figure 1
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$h4("Interactive 3D t-SNE Viewer"),
                          layout_columns(
                            col_widths = c(3, 9),
                            div(
                              selectInput("xcol", "X Axis", choices = c("X", "Y", "Z"), selected = "X"),
                              selectInput("ycol", "Y Axis", choices = c("X", "Y", "Z"), selected = "Y"),
                              selectInput("zcol", "Z Axis", choices = c("X", "Y", "Z"), selected = "Z"),
                              selectInput("color_by", "Color By", choices = c("Cluster"), selected = "Cluster")
                            ),
                            plotOutput("tsne3d", height = "600px")
                          ),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px;",
                                 "Figure 1. 10 clusters determined after analysis on the gene count matrix cleaned with excluding rows with less than 100 counts based on pipeline provided in paper")
                 )
            )
  ),
###Alternative analysis section###  
  nav_panel("Alternative Analysis",
            card("Modern Clustering Pipeline",
                 
#Method section of Alternative analysis
                 tags$div(
                   tags$h4("Methods and Materials"),
                   tags$ul(
                     tags$li(HTML('Trimming reads with <a href="https://github.com/OpenGene/fastp" target="_blank" style="color:#00AEEF;">fastp</a> (v.0.24.0)')),
                     tags$li(HTML('Removing Nextera adapters with fastp</a>')),
                     tags$li(HTML('Indexing transcriptome with <a href="https://github.com/COMBINE-lab/salmon" target="_blank" style="color:#00AEEF;">Salmon</a> (v.1.10.0)')),
                     tags$li(HTML('Quantifying reads with Salmon</a>')),
                     tags$li(HTML('Importing data with <a href="https://bioconductor.org/packages/release/bioc/html/tximport.html" target="_blank" style="color:#00AEEF;">tximport</a> (v.1.28.0)')),
                     tags$li("Creating gene count matrix"),
                     tags$li(HTML('Analyzing matrix with <a href="https://satijalab.org/seurat/" target="_blank" style="color:#00AEEF;">Seurat</a> (v.5.2.1)')),
                     tags$li("Clustering cells (resolution 1.2)"),
                     tags$li("Identifying top 20 expressed genes"),
                     tags$li("Analyzing MHC-I gene expression")
                   )
                 ),
#Challenges of this section
                 tags$div(
                   style = "margin-top: 40px; ",
                   tags$h4("Challenges"),
                   tags$ul(
                     tags$li("Difficulty of merging quantifications obtained from Salmon and converting gene_IDs to gene_names."),
                     tags$li("Finding optimal resolution for clusters in Seurat. "),
                     tags$li("Addressing the duplication in gene names (multiple rows with the same gene name) after loading the gene count matrix"),
                   )
                 ),

#Results of this section
                 tags$div(
                   style = "margin-top: 40px; ",
                   tags$h4("Results"),
                   tags$ul(
                     tags$li("We identified different brain cell clusters representing major cell types."),
                     tags$li("Top 20 most expressed genes were identified in each cluster."),
                     tags$li(HTML('Manual annotation of clusters was carried out, and the results were validated using <a href="https://panglaodb.se/markers.html?cell_type=%27all_cells%27" target="_blank" style="color:#00AEEF;">PanglaoDB</a>.')),
                     tags$li("MHC-I pathway gene expression profiles were investigated in various cells in the dataset.")
                   )
                 ),
                 
                 # Figure 2
                 tags$div(style = "margin-top: 60px; text-align: center;",
                          tags$h4("Figure 2. 3D plot for 9 clusters obtained using alternative analysis"),
                          tags$div(style = "margin: auto; width: 300px;",
                                   sliderInput("angle_alt", "View Angle:", min = 0, max = 180, value = 55)),
                          tags$div(style = "margin: auto; width: 90%;",
                                   plotOutput("scatterplot3d_alt", height = "500px")),
                          tags$div(style = "margin-top: 20px;",
                                   tags$img(src = "Picture3.jpg", style = "max-width: 50%; border-radius: 15px;"))
                 ),
                 
                 
                 # Table 3
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$h5("Table 3. Top 20 most expressed genes in each cluster for the dataset analyzed with Seurat."),
                          div(
                            style = "margin: auto; width: 90%; background-color: #2C2C3A; padding: 20px; border-radius: 15px;",
                            dataTableOutput("top_markers_alt_table")
                          ),
                 ),
                 
                 # Remaining figures and tables
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "file_2025-04-26_18.07.37.png", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 3. 9 clusters found for the dataset analyzed with Seurat.")
                 ),
                 
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "file_2025-04-26_18.14.04.png", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 4. Manual annotation results for pipeline produced with fastp and Seurat.")
                 ),
                 
                 
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "22222222.jpg", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 5. Top expressed genes in each cluster. AQP4: astrocytes, TMEM144: oligodendrocytes, AIF1: microglia, TNR: OPC, SATB2: newly born neurons, TAC1: neuroendocrine cells, CCK: GABAergic neurons, B2M: endothelial cells, TOP2A: replicating neuronal progenitors.")
                 ),
                 
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "file_2025-04-26_18.18.33.png", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 6. S and G2M scores for different clusters. Higher scores indicate that cells in that cluster are more likely to be in replication stage.")
                 ),
                 

                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "V.jpg", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 7. Cluster annotation with PanglaoDB cell type gene expression markers.")
                 ),
                 
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$h5("Table 4. Summary of clusters found with manual annotation and validated with PanglaoDB."),
                          dataTableOutput("last_table"),
                          div(
                            style = "margin: auto; width: 90%; background-color: #2C2C3A; padding: 20px; border-radius: 15px;",
                            dataTableOutput("last_table")
                          ),
                 ),
                 tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "FetalvsAdult.jpg", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 8. Difference between the number of adult and fetal cells obtained from alternative pipeline.")
                 ),
                tags$div(style = "margin-top: 40px; text-align: center;",
                          tags$img(src = "11111.jpg", style = "max-width: 70%; border-radius: 15px;"),
                          tags$p(style = "color: #AAAAAA; margin-top: 10px; font-size: 0.9em; font-style: italic;",
                                 "Figure 9. Expression of different genes involved with MHCI across different cells in the dataset analyzed with alternative tools.")
                 ),
                 
                 
            )
  ),
#Creates new panel for Contact page  
  nav_panel("Contact",
            card("Meet the Team",
                 tags$ul(
                   tags$li("Abirajh Arulrajah - aa1152@student.le.ac.uk"),
                   tags$li("Hesam Moazzen - hm435@student.le.ac.uk"),
                   tags$li("Sanskriti Sanskriti - s82@student.le.ac.uk"),
                   tags$li("Priyadharshini  - pv97@student.le.ac.uk"),
                   tags$li("Alice E. Haskell - aeh35@student.le.ac.uk")
                 )
            )
  )
)

### SERVER section####
server <- function(input, output, session) {
#Logic for the navigation buttons  
  observeEvent(input$go_original, {
    updateNavbarPage(session, inputId = "navbar", selected = "Original Analysis")
  })
  
  observeEvent(input$go_alternative, {
    updateNavbarPage(session, inputId = "navbar", selected = "Alternative Analysis")
  })
  
  observeEvent(input$go_contact, {
    updateNavbarPage(session, inputId = "navbar", selected = "Contact")
  })

#Renders the first table and adds features such as a search bar and highlights 
  output$number_cells_table <- renderDataTable({
    datatable(
      number_cells,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'ftip',
        stripeClasses = c('table-primary', 'table-secondary')
      ),
      class = "display nowrap compact",
      extensions = 'KeyTable'
    )
  })
  
#Renders the second table 
  output$top20_old_table <- renderDataTable({
    datatable(
      top20_old,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'ftip',
        stripeClasses = c('table-primary', 'table-secondary')
      ),
      class = "display nowrap compact",
      extensions = 'KeyTable'
    )
  })
  
#Renders the first plot 
  output$tsne3d <- renderPlot({
    req(input$xcol, input$ycol, input$zcol, input$color_by)
    cluster_levels <- sort(unique(tsne_data[[input$color_by]]))
    colors_palette <- rainbow(length(cluster_levels))
    cluster_colors <- colors_palette[as.numeric(factor(tsne_data[[input$color_by]]))]
    
    scatterplot3d(
      x = tsne_data[[input$xcol]],
      y = tsne_data[[input$ycol]],
      z = tsne_data[[input$zcol]],
      color = cluster_colors,
      pch = 20,
      main = "3D t-SNE Viewer",
      xlab = input$xcol,
      ylab = input$ycol,
      zlab = input$zcol,
      angle = 55
    )
    
    legend("topright", legend = paste("Cluster", cluster_levels),
           col = colors_palette, pch = 20, cex = 1.5, text.col = "black", box.lwd = 0)
  })

#Renders third table  
  output$top_markers_alt_table <- renderDataTable({
    datatable(
      top_markers_alt,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'ftip',
        stripeClasses = c('table-primary', 'table-secondary')
      ),
      class = "display nowrap compact",
      extensions = 'KeyTable'
    )
  })
  
#Renders second plot  
  output$scatterplot3d_alt <- renderPlot({
    scatterplot3d(
      x = scatterplot_data$UMAP_1,
      y = scatterplot_data$UMAP_2,
      z = scatterplot_data$UMAP_3,
      color = as.character(scatterplot_data$Color),
      pch = 20,
      main = "Alternative 3D Scatterplot",
      xlab = "UMAP 1",
      ylab = "UMAP 2",
      zlab = "UMAP 3",
      angle = input$angle_alt
    )
  })

#Renders fourth table  
  output$last_table <- renderDataTable({
    datatable(
      last_table,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'ftip',
        stripeClasses = c('table-primary', 'table-secondary')
      ),
      class = "display nowrap compact",
      extensions = 'KeyTable'
    )
  })
  
}

# Runs App
shinyApp(ui = ui, server = server)
