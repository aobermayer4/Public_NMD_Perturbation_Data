
library(data.table)
library(stringr)
library(GEOquery)
library(dplyr)

setwd("~/R/Shiny/Shiny_App_Paths/")



app_keywords <- fread("App_Path_Parts_Freq_GSE_codes_20260512_v2.txt", na.strings = c("","NA"), header = T, data.table = F)

app_paths <- fread("App_Path_Parts_20260513.txt", na.strings = c("","NA"), header = T, data.table = F)


app_keywords$Source <- gsub("Clincal","Clinical",app_keywords$Source)
app_keywords$Source <- gsub("Rea World","Real World",app_keywords$Source)


app_keywords[which(app_keywords$Project == "Orien"),"Project"] <- "ORIEN"

app_keywords[which(app_keywords$Region == "Hematologic Cancer"),"Region"] <- NA
app_keywords[which(app_keywords$Region == "Prostate Cancer"),"Region"] <- "Prostate"
app_keywords[which(app_keywords$Region == "Colon Cancer"),"Region"] <- "Colon"

app_keywords[which(app_keywords$Disease_Type_Broad == "Colon Cancer"),"Disease_Type_Broad"] <- "Colon/Rectal Cancer"
app_keywords[which(app_keywords$Disease_Type_Broad == "Breast cancer"),"Disease_Type_Broad"] <- "Breast Cancer"

app_keywords[which(app_keywords$Tissue_codes == "Bain"),"Tissue_codes"] <- "Brain"


region_list <- split(app_keywords$Keyword_real,app_keywords$Region)
site_list <- split(app_keywords$Keyword_real,app_keywords$Site)
project_list <- split(app_keywords$Keyword_real,app_keywords$Project)
source_list <- split(app_keywords$Keyword_real,app_keywords$Source)
DisTypeS_list <- split(app_keywords$Keyword_real,app_keywords$Disease_Type_Specific)
DisTypeB_list <- split(app_keywords$Keyword_real,app_keywords$Disease_Type_Broad)
Tissue_list <- split(app_keywords$Keyword_real,app_keywords$Tissue_codes)


anno_list <- list(Project = project_list,
                  Source = source_list,
                  Disease_Region = region_list,
                  Disease_Site = site_list,
                  Disease_Type_Broad = DisTypeB_list,
                  Disease_Type_Specific = DisTypeS_list,
                  Tissue = Tissue_list)

app_paths_anno <- annotate_table(df = app_paths,
                                 string_col = "AppPath",
                                 annotations = anno_list)



write.table(app_paths_anno,"App_Path_Parts_Anno_20260513.txt", sep = '\t', row.names = F)


app_paths_anno2 <- annotate_table(df = app_paths,
                                  string_col = "AppPath",
                                  annotations = anno_list)

colnames(app_paths_anno2)[c(6:11)] <- c("DirLevel1","DirLevel2","DirLevel3","DirLevel4","DirLevel5","DirLevel6")

write.table(app_paths_anno2,"App_Path_Parts_Anno_20260513_v2.txt", sep = '\t', row.names = F)




app_paths_anno3 <- fread("App_Path_Parts_Anno_20260513_v2_working.txt", na.strings = c("","NA"), header = T, data.table = F)


app_paths_anno3_rem <- app_paths_anno3[which(app_paths_anno3$Remove == "x"),]
app_paths_anno4 <- app_paths_anno3[which(is.na(app_paths_anno3$Remove)),]


app_paths_anno5 <- unique(app_paths_anno4[,c(1:14,18,22,29,37,45)])


app_paths_anno6 <- app_paths_anno5[,-c(2,7:12)]



app_paths_anno6$AppLink <- paste0("https://biostools.moffitt.org/4472414/Shiny/",gsub("/share/dept_bbsr/Projects/Shaw_Timothy/Shiny_App_Folder/","",app_paths_anno5$AppPath))


app_paths_anno6 <- app_paths_anno6 %>%
  relocate(AppLink,AppPath,Last_Mod_Date,Last_Mod_Time, .after = last_col()) %>%
  as.data.frame()

app_paths_anno6[which(app_paths_anno6$AppBase == "EASY"),"AppBase"] <- "EASY App"
app_paths_anno6[which(app_paths_anno6$AppBase == "PATH_SURVEYOR"),"AppBase"] <- "PATH SURVEYOR"
app_paths_anno6[which(app_paths_anno6$AppBase == "PATH_SURVEYOR_Lite"),"AppBase"] <- "PATH SURVEYOR Lite"





write.table(app_paths_anno6,"App_Path_Parts_Anno_linked_20260513.txt", sep = '\t', row.names = F)




































