#Fjern alle tex
protect <- c("books","static")
formats <- c(".log",".tex",".toc",".aux")
for (f in formats) {
  files <- list.files(path = "joakim-bilyk.github.io/docs",include.dirs = F, full.names = T, recursive = T,pattern=f)
  for (p in protect) {
    files <- files[!grepl(p,files)]
  }
  file.remove(files)
}
#Fjern alle pdf i root
files <- list.files(path = "joakim-bilyk.github.io/docs",include.dirs = F, full.names = T, recursive = F,pattern=".pdf")
for (p in protect) {
  files <- files[!grepl(p,files)]
}
file.remove(files)
files <- list.files(path = "joakim-bilyk.github.io/docs/pdf/render",include.dirs = T, full.names = T)
files <- files[grepl("_files",files)]
unlink(files,recursive = TRUE)

formats <- c(".log",".tex")
for (f in formats) {
  files <- list.files(include.dirs = F, full.names = T, recursive = F,pattern=f)
  files <- files[!grepl("my-template.tex",files)]
  for (p in protect) {
    files <- files[!grepl(p,files)]
  }
  file.remove(files)
}

# Render theory book
files <- list.files(path = "_assets/theory/docs",include.dirs = T, full.names = T)
file.remove(files)
rmarkdown::render_site("_assets/theory",encoding = 'ISO8859-1')
unlink("joakim-bilyk.github.io/docs/books",recursive = TRUE)
dir.create("joakim-bilyk.github.io/docs/books")
dir.create("joakim-bilyk.github.io/docs/books/theory")
files <- list.files(path = "_assets/theory/docs",include.dirs = T, full.names = T)
file.copy(files,"joakim-bilyk.github.io/docs/books/theory",recursive = TRUE)
file.remove("joakim-bilyk.github.io/docs/pdf/render/theory.pdf")
file.copy("joakim-bilyk.github.io/docs/books/theory/theory.pdf","joakim-bilyk.github.io/docs/pdf/render/theory.pdf")

#Render exercise book
files <- list.files(path = "_assets/exercises/docs",include.dirs = T, full.names = T)
file.remove(files)
rmarkdown::render_site("_assets/exercises",encoding = 'ISO8859-1')
dir.create("joakim-bilyk.github.io/docs/books/exercises")
files <- list.files(path = "_assets/exercises/docs",include.dirs = T, full.names = T)
file.copy(files,"joakim-bilyk.github.io/docs/books/exercises",recursive = TRUE)
file.remove("joakim-bilyk.github.io/docs/pdf/render/exercises.pdf")
file.copy("joakim-bilyk.github.io/docs/books/exercises/exercises.pdf","joakim-bilyk.github.io/docs/pdf/render/exercises.pdf")

#Render exam prep book
files <- list.files(path = "_assets/exam/docs",include.dirs = T, full.names = T)
file.remove(files)
rmarkdown::render_site("_assets/exam",encoding = 'ISO8859-1')
dir.create("joakim-bilyk.github.io/docs/books/exam")
files <- list.files(path = "_assets/exam/docs",include.dirs = T, full.names = T)
file.copy(files,"joakim-bilyk.github.io/docs/books/exam",recursive = TRUE)
file.remove("joakim-bilyk.github.io/docs/pdf/render/exam.pdf")
file.copy("joakim-bilyk.github.io/docs/books/exam/exam.pdf","joakim-bilyk.github.io/docs/pdf/render/exam.pdf")

#Get statics
unlink("joakim-bilyk.github.io/docs/static",recursive = TRUE)
dir.create("joakim-bilyk.github.io/docs/static")
files <- list.files(path = "_assets/static",include.dirs = T, full.names = T)
file.copy(files,"joakim-bilyk.github.io/docs/static",recursive = TRUE)
