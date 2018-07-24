
## READS:
## dumped nematostella_targets from 2 mysqldev instances 10 times each
## e.g.: "time mysqldump -h mysqldev nematostella_targets > nt1.sql", using 'real' time in seconds


bytes.old <- 198371656
bytes.new <- 198371559

bytes.old/1024^2


new <- c(5.468,5.508,5.474,5.376,5.464,5.482,5.511,5.568,5.461,5.415)  # mysqldev
old <- c(6.513,6.214,6.264,5.955,6.109,6.398,6.450,6.271,6.418,6.378)  # mysql-dev

mean(new); sd(new); (bytes.new/1024^2)/mean(new)
mean(old); sd(old); (bytes.old/1024^2)/mean(old)

((bytes.new/1024^2)/mean(new)) / ((bytes.old/1024^2)/mean(old))





## WRITES:
## pushed nematostella_targets dumps to mysqldev instances 11 times each.
## the first push in each test set was ignored.  Databases were dropped and re-created between each set.
## e.g.: "time mysql -h mysql-dev apa_nemato_test < nt-1.sql", using 'real' time in seconds.
## mysqldev received both mysql-dev and mysqldev dumps (new.old vs new.new).
## mysql-dev received only mysql-dev dumps (old.old).


new.old <- c(14.491,14.362,14.502,14.420,14.502,14.501,14.388,14.326,14.270,14.211)
new.new <- c(14.157,14.246,14.215,14.323,14.250,14.208,14.183,14.287,14.293,14.247)
old.old <- c(27.857,27.553,28.127,28.404,27.876,27.687,28.089,26.728,27.793,27.622)

mean(new.new); sd(new.new); (bytes.new/1024^2)/mean(new.new)
mean(new.old); sd(new.old); (bytes.new/1024^2)/mean(new.old)
mean(old.old); sd(old.old); (bytes.old/1024^2)/mean(old.old)

((bytes.new/1024^2)/mean(new.new)) / ((bytes.old/1024^2)/mean(old.old))
((bytes.new/1024^2)/mean(new.old)) / ((bytes.old/1024^2)/mean(old.old))

