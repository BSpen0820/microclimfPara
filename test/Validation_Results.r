
library(tools)
library(digest)

files <- list.files(path = "./", pattern = "\\.rds$", full.names = TRUE, recursive = T) 

files.df <- data.frame(
  file = files,
  name = file_path_sans_ext(basename(files)),
  dir = dirname(files)
)

files.df$hash <- NA

for(i in 1:nrow(files.df)) {
  files.df$hash[i] <- digest(files.df$file[i], algo = "md5", file = TRUE)
}

files.df$identical_2_Base <- NA
files.df$identical_Toler_Base <- NA

for(i in 1:5){
    #i <- 3
    unique_name <- unique(files.df$name)[i]

    files_subset <- files.df[files.df$name == unique_name, ]
    base_hash <- files_subset$hash[files_subset$dir == "./Base_Results"]

    for(j in 1:3){
        #j <- 2
        if(files_subset$dir[j] == "./Base_Results") next
        compare_hash <- files_subset$hash[files_subset$dir == files_subset$dir[j]]
        same <- base_hash == compare_hash
        files.df$identical_2_Base[files.df$file == files_subset$file[j]] <- same
        if(same) {
            files.df$identical_Toler_Base[files.df$file == files_subset$file[j]] <- TRUE
        } else {
            base_data <- readRDS(files_subset$file[files_subset$dir == "./Base_Results"])
            compare_data <- readRDS(files_subset$file[j])
            # learned that cpp with parallel worker will return values a little different due to passing values with "less" precision,
            # but essentially the same results, so I set a tolerance level to check for "near" equality.
            files.df$identical_Toler_Base[files.df$file == files_subset$file[j]] <- all.equal(base_data, compare_data, tolerance = 1e-10) 
        }

    }

}

files.df
write.csv(files.df, "Validation_Results.csv", row.names = FALSE)
