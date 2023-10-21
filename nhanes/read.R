library(dplyr)
library(readr)

pa = "/home/kshedden/data/Teaching/nhanes/2017-2018"

fn = c("DEMO_J.csv.gz", "BMX_J.csv.gz", "BPX_J.csv.gz")

df = NULL
for (f in fn) {
    dx = read_csv(sprintf("%s/%s", pa, f))
    if (is.null(df)) {
        df = dx
    } else {
        df = left_join(df, dx, "SEQN")
    }
}

df = df %>% mutate(RIAGENDR=recode(df$RIAGENDR, `1`="M", `2`="F"))
df = df %>% mutate(RIDRETH1=recode(df$RIDRETH1, `1`="MA", `2`="OH", `3`="NHW", `4`="NHB", `5`="Other"))

df = df %>% filter(RIDAGEYR >= 18)
