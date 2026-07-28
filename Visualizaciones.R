library(readr)
library(ggpubr)

#### Graficos Jaccard de PCA + k-means####
# Algunos gráficos descriptivos de los indices de Jaccard generados en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.

jac12 <- read_table2("Objetos/Jaccard_ind_12.csv")[,2] 
jac22 <- read_table2("Objetos/Jaccard_ind_22.csv")[,2]
jac32 <- read_table2("Objetos/Jaccard_ind_32.csv")[,2]
jac11 <- read_table2("Objetos/Jaccard_ind_11.csv")[,2] 
jac21 <- read_table2("Objetos/Jaccard_ind_21.csv")[,2] 
jac31 <- read_table2("Objetos/Jaccard_ind_31.csv")[,2]

jacc1 <- as.data.frame(jac12)
jacc2 <- as.data.frame(jac22)
jacc3 <- as.data.frame(jac32)
mean(jacc1$`"x"`)
mean(jacc2$`"x"`)
mean(jacc3$`"x"`)
#cbind(type=1, jac11),cbind(type=4, jac31),
df <- rbind(cbind(type=11, jac11),cbind(type=12, jac12),
            cbind(type=31, jac31),cbind(type=32, jac32),
            cbind(type=22, jac22),cbind(type=21, jac21)
            )
names(df) <- c("type", "index")
head(df)

## Boxplot del indice de Jaccard
df %>%
  ggboxplot(x       = "type", 
            y       = "index",
            fill    = "type",
            palette = c("#27763D","#86CA78",
                        "#C71F3C","#F8806A", 
                        "#32769B","#63B1C1"))

## Violin del indice de Jaccard
df %>%
  ggviolin(x          = "type", 
           y          = "index", 
           xlab       = "Tipo de simulación",
           ylab       = "Puntaje de Jaccard",
           fill       = "type",
           alpha      = 0.8,
           palette    = c("#27763D","#86CA78",
                          "#C71F3C","#F8806A", 
                          "#32769B","#63B1C1"),
           add        = "boxplot")

## Densidad
df$type <- as.factor(df$type)
df %>%
  ggdensity(x       = "index", 
            add     = "mean", 
            rug     = TRUE,
            color   = "type", 
            fill    = "type",
            xlab    = "Puntaje de Jaccard",
            palette =  c("#27763D","#86CA78",
                         "#C71F3C","#F8806A", 
                         "#32769B","#63B1C1"))

## Histograma
df %>%
  gghistogram(x       = "index",
              add     = "mean", 
              rug     = TRUE,
              color   = "type", 
              fill    = "type",
              xlab    = "Puntaje de Jaccard",
              palette = c("#27763D","#86CA78",
                          "#C71F3C","#F8806A", 
                          "#32769B","#63B1C1"))

## Densidad e histograma
p1 <- df %>%
  ggdensity(x       = "index", 
            add     = "mean", 
            rug     = TRUE,
            color   = "type", 
            fill    = "type",
            xlab    = "Puntaje de Jaccard",
            palette =  c("#2e00fa","#a000bc","#ca0086")
  )
p2 <- df %>%
  gghistogram(x       = "index",
              add     = "mean", 
              rug     = TRUE,
              color   = "type", 
              fill    = "type",
              xlab    = "Puntaje de Jaccard",
              palette = c("#2e00fa","#a000bc","#ca0086")
  )
ggarrange(p1, 
          p2, 
          labels = c("A", "B"),
          ncol   = 2
)

#### Graficos Varianza retenida de PCA + k-means####
# Algunos gráficos descriptivos de la varianza retenida generadas en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.

var1 <- read_table2("Objetos/Varianza_ret_12.csv")[,2]
var2 <- read_table2("Objetos/Varianza_ret_22.csv")[,2]
var3 <- read_table2("Objetos/Varianza_ret_32.csv")[,2]
var4 <- read_table2("Objetos/Varianza_ret_14.csv")[,2]
var5 <- read_table2("Objetos/Varianza_ret_24.csv")[,2]
var6 <- read_table2("Objetos/Varianza_ret_34.csv")[,2]
df <- rbind(cbind(type=1, var1),cbind(type=2, var2),cbind(type=3, var3), cbind(type=1, var4),cbind(type=2, var5),cbind(type=3, var6))
names(df) <- c("type", "var")
head(df)

## Boxplot del indice de Jaccard
df %>%
  ggboxplot(x       = "type", 
            y       = "var",
            fill    = "type",
            palette = c("#2e00fa", "#a000bc", "#ca0086"))

## Violin del indice de Jaccard
df %>%
  ggviolin(x          = "type", 
           y          = "var", 
           xlab       = "Tipo de simulación",
           ylab       = "Varianza retenida",
           fill       = "type",
           alpha      = 0.8,
           palette    = c("#2e00fa","#a000bc","#ca0086"),
           add        = "boxplot")

## Densidad
df$type <- as.factor(df$type)
df %>%
  ggdensity(x       = "var", 
            add     = "mean", 
            rug     = TRUE,
            color   = "type", 
            fill    = "type",
            xlab    = "Varianza retenida",
            palette =  c("#2e00fa","#a000bc","#ca0086"))

## Histograma
df %>%
  gghistogram(x       = "var",
              add     = "mean", 
              rug     = TRUE,
              color   = "type", 
              fill    = "type",
              xlab    = "Varianza retenida",
              palette = c("#2e00fa","#a000bc","#ca0086"))

## Densidad e histograma
p1 <- df %>%
  ggdensity(x       = "var", 
            add     = "mean", 
            rug     = TRUE,
            color   = "type", 
            fill    = "type",
            xlab    = "Varianza Retenida",
            palette =  c("#2e00fa","#a000bc","#ca0086")
  )
p2 <- df %>%
  gghistogram(x       = "var",
              add     = "mean", 
              rug     = TRUE,
              color   = "type", 
              fill    = "type",
              xlab    = "Varianza Retenida",
              palette = c("#2e00fa","#a000bc","#ca0086")
  )
ggarrange(p1, 
          p2, 
          labels = c("A", "B"),
          ncol   = 2
)

#### Graficos Componentes retenidas de PCA + k-means####
# Algunos gráficos descriptivos de las componentes retenidas generadas 
# en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.

comp1 <- read_table2("Objetos/Componentes_ret_12.csv")[,2]
comp2 <- read_table2("Objetos/Componentes_ret_22.csv")[,2]
comp3 <- read_table2("Objetos/Componentes_ret_32.csv")[,2]
df <- rbind(cbind(type=1, comp1),cbind(type=2, comp2),cbind(type=3, comp3))
names(df) <- c("type", "comp")
head(df)

## Boxplot del indice de Jaccard
df %>%
  ggboxplot(x       = "type", 
            y       = "comp",
            fill    = "type",
            palette = c("#2e00fa", "#a000bc", "#ca0086"))

## Violin del indice de Jaccard
df %>%
  ggviolin(x          = "type", 
           y          = "comp", 
           xlab       = "Tipo de simulación",
           ylab       = "Varianza retenida",
           fill       = "type",
           alpha      = 0.8,
           palette    = c("#2e00fa","#a000bc","#ca0086"),
           add        = "boxplot")

#### Graficos Jaccard de UMAP + k-means####
# Algunos gráficos descriptivos de los indices de Jaccard generados en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.
library(dplyr)

jac1 <- read_table2("Objetos/Jaccard_ind_11.csv")[,2] 
jac2 <- read_table2("Objetos/Jaccard_ind_21.csv")[,2]
jac3 <- read_table2("Objetos/Jaccard_ind_31.csv")[,2]
jacc1 <- as.data.frame(jac1)
jacc2 <- as.data.frame(jac2)
jacc3 <- as.data.frame(jac3)
mean(jacc1$`"x"`)
mean(jacc3$`"x"`)
mean(jacc2$`"x"`)

Coords_UMAP_11 <- read_table2("Objetos/Coords_UMAP_11.csv") %>% as.data.frame()
Cluster_generados_11 <- read_table2("Objetos/Cluster_generados_11.csv") %>% as.data.frame()
plot(Coords_UMAP_11[,2:3], col=Cluster_generados_11[,21], xlab="UMAP 1", ylab="UMAP2")
Cluster_originales_11 <- read_table2("Objetos/Cluster_originales_11.csv") %>% as.data.frame()
plot(Coords_UMAP_11[,2:3], col=Cluster_originales_11[,21], xlab="UMAP 1", ylab="UMAP2")

Coords_UMAP_21 <- read_table2("Objetos/Coords_UMAP_21.csv") %>% as.data.frame()
Cluster_generados_21 <- read_table2("Objetos/Cluster_generados_21.csv") %>% as.data.frame()
plot(Coords_UMAP_21[,2:3], col=Cluster_generados_21[,21], xlab="UMAP 1", ylab="UMAP2")
Cluster_originales_21 <- read_table2("Objetos/Cluster_originales_21.csv") %>% as.data.frame()
plot(Coords_UMAP_21[,2:3], col=Cluster_originales_21[,21], xlab="UMAP 1", ylab="UMAP2")

Coords_UMAP_31 <- read_table2("Objetos/Coords_UMAP_31.csv") %>% as.data.frame()
Cluster_generados_31 <- read_table2("Objetos/Cluster_generados_31.csv") %>% as.data.frame()
plot(Coords_UMAP_31[,2:3], col=Cluster_generados_31[,21], xlab="UMAP 1", ylab="UMAP2")
plot(Coords_UMAP_31[,2:3], xlab="UMAP 1", ylab="UMAP2")
Cluster_originales_31 <- read_table2("Objetos/Cluster_originales_31.csv") %>% as.data.frame()
plot(Coords_UMAP_31[,2:3], col=Cluster_originales_31[,21], xlab="UMAP 1", ylab="UMAP2")

#### Graficos hyper de UMAP####
# Algunos gráficos descriptivos de la varianza retenida generadas en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.

vec1 <- read_table2("Objetos/Hyper_UMAP_11.csv")[2,-1] 
vec2 <- read_table2("Objetos/Hyper_UMAP_21.csv")[2,-1]
vec3 <- read_table2("Objetos/Hyper_UMAP_31.csv")[2,-1]
vec4 <- read_table2("Objetos/Hyper_UMAP_13.csv")[2,-1]
vec5 <- read_table2("Objetos/Hyper_UMAP_23.csv")[2,-1]
vec6 <- read_table2("Objetos/Hyper_UMAP_33.csv")[2,-1]
c(t(vec1), t(vec4))
t(vec1)
df <- rbind(cbind(type=1, t(vec1)), cbind(type=1, t(vec4)), cbind(type=2, t(vec2)), cbind(type=2, t(vec5)), cbind(type=3, t(vec3)), cbind(type=3, t(vec6)))
df <- as.data.frame(df)
names(df) <- c("tipo", "vecinos")
head(df)

## Boxplot del indice de Jaccard
df %>%
  ggboxplot(x       = "tipo", 
            y       = "vecinos",
            fill    = "tipo",
            palette = c("#2e00fa", "#a000bc", "#ca0086"))

dis1 <- read_table2("Objetos/Hyper_UMAP_11.csv")[3,-1] 
dis2 <- read_table2("Objetos/Hyper_UMAP_21.csv")[3,-1]
dis3 <- read_table2("Objetos/Hyper_UMAP_31.csv")[3,-1]
dis4 <- read_table2("Objetos/Hyper_UMAP_13.csv")[3,-1]
dis5 <- read_table2("Objetos/Hyper_UMAP_23.csv")[3,-1]
dis6 <- read_table2("Objetos/Hyper_UMAP_33.csv")[3,-1]
c(t(dis1), t(dis4))
t(dis1)
df <- rbind(cbind(type=1, t(dis1)), cbind(type=1, t(dis4)), cbind(type=2, t(dis2)), cbind(type=2, t(dis5)), cbind(type=3, t(dis3)), cbind(type=3, t(dis6)))
df <- as.data.frame(df)
names(df) <- c("tipo", "distancia")
head(df)

## Boxplot del indice de Jaccard
df %>%
  ggboxplot(x       = "tipo", 
            y       = "distancia",
            fill    = "tipo",
            palette = c("#2e00fa", "#a000bc", "#ca0086"))


#### Graficos Jaccard de PCA + Louvain####
# Algunos gráficos descriptivos de los indices de Jaccard generados en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.

jac14 <- read_table2("Objetos/Jaccard_ind_14.csv")[,2] 
jac24 <- read_table2("Objetos/Jaccard_ind_24.csv")[,2]
jac34 <- read_table2("Objetos/Jaccard_ind_34.csv")[,2]
jac13 <- read_table2("Objetos/Jaccard_ind_13.csv")[,2] 
jac23 <- read_table2("Objetos/Jaccard_ind_23.csv")[,2] 
jac33 <- read_table2("Objetos/Jaccard_ind_33.csv")[,2]

jacc14 <- as.data.frame(jac14)
jacc24 <- as.data.frame(jac24)
jacc34 <- as.data.frame(jac34)
mean(jacc14$`"x"`)
mean(jacc24$`"x"`)
mean(jacc34$`"x"`)
#cbind(type=1, jac11),cbind(type=4, jac31),

df <- rbind(cbind(type=13, jac13),cbind(type=14, jac14),
            cbind(type=33, jac33),cbind(type=34, jac34),
            cbind(type=23, jac23),cbind(type=24, jac24)
)
names(df) <- c("type", "index")
head(df)

## Boxplot del indice de Jaccard
df %>%
  ggboxplot(x       = "type", 
            y       = "index",
            fill    = "type",
            palette = c("#27763D","#86CA78",
                        "#C71F3C","#F8806A", 
                        "#32769B","#63B1C1"))

## Violin del indice de Jaccard
df %>%
  ggviolin(x          = "type", 
           y          = "index", 
           xlab       = "Tipo de simulación",
           ylab       = "Puntaje de Jaccard",
           fill       = "type",
           alpha      = 0.8,
           palette    = c("#27763D","#86CA78",
                          "#C71F3C","#F8806A", 
                          "#32769B","#63B1C1"),
           add        = "boxplot")

## Densidad
df$type <- as.factor(df$type)
df %>%
  ggdensity(x       = "index", 
            add     = "mean", 
            rug     = TRUE,
            color   = "type", 
            fill    = "type",
            xlab    = "Puntaje de Jaccard",
            palette =  c("#27763D","#86CA78",
                         "#C71F3C","#F8806A", 
                         "#32769B","#63B1C1"))

## Histograma
df %>%
  gghistogram(x       = "index",
              add     = "mean", 
              rug     = TRUE,
              color   = "type", 
              fill    = "type",
              xlab    = "Puntaje de Jaccard",
              palette = c("#27763D","#86CA78",
                          "#C71F3C","#F8806A", 
                          "#32769B","#63B1C1"))

#### Graficos Jaccard de UMAP + Louvain####
# Algunos gráficos descriptivos de los indices de Jaccard generados en las tres
# simulaciones, reducidos con PCA y agrupados con k-means.
library(dplyr)

jac13 <- read_table2("Objetos/Jaccard_ind_13.csv")[,2] 
jac23 <- read_table2("Objetos/Jaccard_ind_23.csv")[,2]
jac33 <- read_table2("Objetos/Jaccard_ind_33.csv")[,2]
jacc13 <- as.data.frame(jac13)
jacc2 <- as.data.frame(jac23)
jacc3 <- as.data.frame(jac33)
mean(jacc13$`"x"`)
mean(jacc3$`"x"`)
mean(jacc2$`"x"`)
