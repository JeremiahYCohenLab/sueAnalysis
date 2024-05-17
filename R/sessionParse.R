session = "689514_2024-02-01_18-06-43"
root = "F:/acuteBehavior"
split_string = strsplit(session, "_")
aniID = split_string[[1]][1]
date = split_string[[1]][2]
time = split_string[[1]][3]

path = gsub("\\\\", "/", root)

# Print the modified path
print(path)

sessionPath = file.path(root, aniID, session)

sortedFolder = file.path(sessionPath, "sorted")


# Directory path

# Find files matching the pattern
matching_files <- list.files(sortedFolder, pattern = paste0("^", session, ".*", '.nwb'), full.names = TRUE)

# Print matching file names
print(matching_files)


