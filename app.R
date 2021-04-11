library(shiny)
library(tidyverse)
library(devtools)
library(rawDiag) #https://github.com/fgcz/rawDiag/releases
library(protViz)
library(plotly)
library(lubridate)
library(scales)
library(data.table)
library(RColorBrewer)
library(markdown)

# Set Shiny options
options(shiny.maxRequestSize=900*1024^2)

# Import amino acid residue mass table
# Masses from http://www.matrixscience.com/help/aa_help.html

masses <- tibble::tribble(
  ~Residue, ~MI_mass,
       "A",  71.0371,
       "C", 103.0092,
       "D", 115.0269,
       "E", 129.0426,
       "F", 147.0684,
       "G",  57.0215,
       "H", 137.0589,
       "I", 113.0841,
       "K",  128.095,
       "L", 113.0841,
       "M", 131.0405,
       "N", 114.0429,
       "P",  97.0528,
       "Q", 128.0586,
       "R", 156.1011,
       "S",   87.032,
       "T", 101.0477,
       "V",  99.0684,
       "W", 186.0793,
       "Y", 163.0633
  )

fl <- list.files(pattern = "*.md")

# 2. UI -----
ui <- fluidPage(
    theme = shinythemes::shinytheme("cerulean"), #https://shiny.rstudio.com/gallery/shiny-theme-selector.html
    
    # Application title
    titlePanel("Mass Spec - Targeted Analysis"),
    headerPanel(""),
    
    #Sidebar (INPUTS)
    sidebarLayout(
        sidebarPanel(
            textAreaInput("targets", "Peptide Targets", ""),
            strong("Raw files"),
            br(),br(),
            tags$div(class="input-group", style="height: 38px;",
                     list(
                         tags$label(class="input-group-btn", style="height: 38px;",
                                    list(
                                        tags$span(class="btn btn-default btn-file", style="top: 0px; bottom: 0px;","Browse...",
                                                  list(
                                                      tags$input(id="csvs", name="csvs", type="file", style="display: none;", multiple="multiple", class="shiny-bound-input")))),
                                    tags$input(type="text", class="form-control", placeholder="No file selected", readonly="readonly",
                                               style="height: 118px; padding-bottom: 8px; top: 0px;")))),
            tags$div(id="csvs_progress", class="progress progress-striped active shiny-file-input-progress",
                     list(
                         tags$div(class="progress-bar"))),
            # Button
            downloadButton("downloadData", "Download")
        ),
        
        
        # Main panel (OUTPUTS)
        mainPanel(
            tabsetPanel(
                tabPanel("README",
                         h2("Popular Targets"),
                         p("You can copy-paste the targets below into the 'Peptide Targets' window to the top left
                         of the screen. To enter your own targets, note the formatting required (mass shifts of PTMs
                           between parentheses and charge state between brackets appended at the end)."),
                         fluidRow(
                           column(4,
                           h3("BSA",style="color:black"),
                           h4("EAC(+57.02)FAVEGPK[+2]",style="color:black"),
                           h4("YLYEIAR[+2]",style="color:black"),
                           h4("DDPHAC(+57.02)YSTVFDK[+2]",style="color:black"),
                           h4("FSALTPDETYVPK[+2]",style="color:black"),
                           h4("LVNELTEFAK[+2]",style="color:black"),
                           h4("QTALVELLK[+2]",style="color:black")),
                           
                           column(4,
                            h3("MC38-Global",style="color:black"),
                            h4("YGYSNRVVDLM[+2]",style="color:black")),
                           
                           column(4,
                            h3("MC38-MHCI",style="color:black"),
                            h4("ASMTNMELM[+2]",style="color:black"),
                            h4("ASMTNM(+15.99)ELM[+2]",style="color:black"),
                            h4("AQLANDVVL[+2]",style="color:black"),
                            h4("SIIVFNLL[+2]",style="color:black"),
                            h4("SSPYSLHYL[+2]",style="color:black"),
                            h4("AALLNSAVL[+2]",style="color:black"),
                            h4("MAPIDHTTM[+2]",style="color:black"),
                            h4("DPSSSVLFEY[+2]",style="color:black"))
                         ),
                         
                        
                         h2("Instructions"),
                         p("In order to get the script to run you need to do two things, in this order: First, indicate
                           which targets you want to look for by copy-pasting text or writing in the text window above. 
                           Second, indicate which Thermo .raw files you want to look for the targets in by dragging and
                           dropping the files into the window (or using the browse button). If you need a sample
                           .raw files, https://www.ebi.ac.uk/pride/archive/projects/PXD009542 is a really nice
                           PRIDE entry.")
                ),
                tabPanel("Targets",
                         tableOutput("table")
                ),
                tabPanel("Retention",
                         plotlyOutput("plot")
                )
            ))))

# 3. Server -----
server <- function(input,output,session) {
    
    # * 3.1 Observe Upload ----- 
    peptide_info <- reactiveValues(info_peptide = NULL)
    file_df <- reactiveValues(df_files = NULL)
    mass_list <- reactiveValues(list_mass = NULL)
    
    observeEvent(input$csvs, {
        info_peptide <- read.table(text = input$targets,header=F,stringsAsFactors = F)
        
        # Breaking up strings into component parts
        broken_list <- list()
        for (a in 1:nrow(info_peptide)){
            broken_list[a] <- strsplit(info_peptide$V1[a],"(?<=.)(?=[\\[])", perl = TRUE)
        }
        
        
        # Now we need to cycle through the broken list one element at a time
        # Each temp vector will represent a sequence (e.g. [[1]])
        # This master loop has a single empty vector to fill, total_mass
        total_mass <- vector()
        
        for (i in 1:length(broken_list)){
            
        }
        # Then, we need to have two separate loops embedded within the one-vector-at at time loop
        # The first loop needs three empty vectors, 
        #temp_mass_mod
        #temp_res_mass
        #temp_combined
        # The first loop goes through each residue in the vector looking for mass mods
        # If it finds a mass mod, it stores that in the temp_mass_mod vector
        # It then removes all mass mods
        # The second loop does a mass lookup on all remaining residues
        # These masses are then combined in the combo vector for each residue
        # Then, outside the second loop but before moving onto the next sequence
        # b-ions and y-ions are calculated, with a fixed number of each b1 to b25
        # and y1 to y25
        # These are assumed to be singly charged
        # Then the total mass is calculated by summing the combo_vector for each sequence
        # When every sequence has been run, the total mass and b1-25 + y1-25 ions are compiled
        # Into a big dataframe
        # m/z is then calculated below like it was for everything else
        
        
        # Extract everything between parentheses out of your vector
        to_add_as_list <- regmatches(info_peptide$V1,gregexpr("(?<=\\().*?(?=\\))", 
                                                              info_peptide$V1, perl=TRUE))
        # Now to_add is a list where to_add[[n]] is the row
        # Coerce to numeric and then sum
        to_add <- vector()
        for (i in 1:length(to_add_as_list)){
            temp <- sum(as.numeric(to_add_as_list[[i]]))
            to_add <- c(to_add,temp)
        }
        
        # Do the same as above, but with charge states at the end e.g. [+2]
        charge_list <- regmatches(info_peptide$V1,gregexpr("(?<=\\[).*?(?=\\])", 
                                                           info_peptide$V1, perl=TRUE))
        
        charge <- vector()
        for (j in 1:length(charge_list)){
            temp <- as.numeric(charge_list[[j]])
            charge <- c(charge,temp)
        }
        
        #Once you've extracted values you can sub everything between parentheses
        clean_names <- info_peptide$V1
        clean_names <- gsub("(?<=\\().*?(?=\\))","",clean_names,perl=T)
        clean_names <- gsub("(?<=\\[).*?(?=\\])","",clean_names,perl=T)
        
        #Sub out the parentheses themselves
        clean_names <- gsub("\\(","",clean_names,perl=T)
        clean_names <- gsub("\\)","",clean_names,perl=T)
        clean_names <- gsub("\\[","",clean_names,perl=T)
        clean_names <- gsub("\\]","",clean_names,perl=T)
        
        # Then, for each element of the vector, have a lookup table with the mass
        # and replace the letter with the residue mass
        
        summed_residues <- vector()
        
        for (n in 1:length(clean_names)){
            temp <- clean_names[n]
            temp <- strsplit(temp,"")
            temp <- temp[[1]]
            for (m in 1:length(temp)){
                temp[m] <- gsub(temp[m],masses[which(masses==temp[m]),2],temp[m])
            }
            temp_sum <- sum(as.numeric(temp))
            summed_residues <- c(summed_residues,temp_sum)
        }
        
        # Then sum the between-parentheses values, residue masses, and add
        # N and C termini and that's it
        total_sum <- vector()
        ppm <- 0 #correction value (0.00047 for my BSA data)
        
        for (j in 1:length(clean_names)){
            total_sum <- c(total_sum,
                           sum(to_add[j],summed_residues[j],18.0106, # Monoisotopic mass of water is 18.01056Da according to https://dodona.ugent.be/en/exercises/575192986/
                               ppm)) #ppm correction value
        }
        
        # Calculate m/z value
        mz <- vector()
        
        for (k in 1:length(total_sum)){
            mz <- c(mz,
                    ((total_sum[k]+(charge[k]*1.0078))/charge[k])
            )
        }
        
        
        
        info_peptide <- data.frame(info_peptide$V1,
                                   to_add,
                                   summed_residues,
                                   total_sum,
                                   charge,
                                   mz)
        
        colnames(info_peptide) <- c("peptide","mods","summed_residues","M","charge","mz")
        
        peptide_info$info_peptide <- info_peptide
        
        
        # ** 3.1.1 File Names -----
        file_df$df_files <- as.data.frame(input$csvs)
        
        # ** 3.1.4 Mass List -----
        list_mass <- data.frame()
        withProgress(message = 'Processing', value = 0, {
            for (i in 1:nrow(file_df$df_files)) {
                incProgress(i/nrow(file_df$df_files), detail = paste0("File ",i," of ",nrow(file_df$df_files)))
                temp_TIC <- read.raw(file_df$df_files[i,4])
                for (j in 1:nrow(info_peptide)) {
                    temp <- readXICs(rawfile=file_df$df_files[i,4],masses=info_peptide$mz[j])
                    df_temp <- data.frame(times=temp[[1]][2],
                                          intensities=temp[[1]]$intensities/(as.numeric(sum(temp_TIC$TIC))/1e11)) %>%
                        mutate(file=file_df$df_files[i,1],peptide=info_peptide$peptide[j])
                    list_mass <- rbind(list_mass,df_temp)
                }
            }
        })
        
        mass_list$list_mass <- list_mass
    })
    
    # * 3.2 Tab1 (Targets) -----  
    output$table <- renderTable({
        peptide_info$info_peptide
        
    },digits=4)
    
    # * 3.3 Tab2 (Retention) -----
    output$plot <- renderPlotly({
        
        plot_ly(mass_list$list_mass,x=~times,y=~intensities, color=~file,
                mode="lines+markers",
                type="scatter",
                text = ~peptide,
                hovertemplate = paste(
                    "<b>%{text}</b><br>",  #<br><br>",
                    "%{xaxis.title.text}: %{x:.2f}<br>",
                    "%{yaxis.title.text}: %{y:.0e}<br>",
                    "<extra></extra>")
        ) %>%
            layout(
                # xaxis = list(
                #   range = c(10,20)),
                yaxis = list(
                    exponentformat = "E")
            )
        
    })
    
    # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      gsub(":","-",gsub(" ","_",paste0(now(),"_targeted_ms.csv")))
    },
    content = function(file) {
      write.csv(mass_list$list_mass, file, row.names = FALSE)
    }
  )
}
app <- shinyApp(ui = ui, server = server)
runApp(app, host ="0.0.0.0", port = 3838, launch.browser = TRUE)

####################
### FUTURE IDEAS ###
####################

# How to work up the MS2 plots
# Store a variable for MS2 scan number as null
# On a click from the user, all MS2 scans within -/+ 3mins that match -/+5ppm
# Of the target peptide are chosen, then within those the highest intensity
# scan is chosen and the MS2 scan number variable is set
# Once that variable is set a labeled graph appears below the RT graph
# Then, for that peptide, its b- and y- ions are called up
# If a mass matches any b- and y-ions (within 2ppm?) above a certain intensity threshold
# (2.5%?) then it is colored green for match and grey for no match
# More complex coloration to indicate confidence of a good match, etc. can be included later easily

# It might also be good to make fixed and variable options and charge states,
# allow the program to make all possible combos and search for them, and 
# then if the output is below a certain value, remove it from the output so
# it doesn't get too crowded

###################
### GARBAGE BIN ###
###################

# d <- strsplit(a, "(?<=.)(?=[A-Z])",perl = TRUE)[[1]]

# length(strsplit(test,"[A-Z]",perl=T)[[1]])
# test2 <- strsplit(test,"\\(",perl=T)[[1]]
# as.numeric(sub("\\)","",test2[2]))