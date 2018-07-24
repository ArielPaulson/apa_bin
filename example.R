#!/usr/bin/env Rscript
 
self.file <- "/home/apa/local/bin/example.R"
self.code <- scan(self.file, what="", sep="\n")
self.code <- c(
    self.code,
    "print(0)",
    "message('\nI just ran an extra command!\n')"
)
 
message("\nI am editing myself...\n")
write.table(self.code, self.file, quote=FALSE, row.names=FALSE, col.names=FALSE)
 
print(1:10)
message("\nThis is EOF\n")
##EOF
