library(seqinr)
library(knitr)


mx_2020 = read.fasta("mx_2020_b519.fasta")
mx_2021 = read.fasta("mx_2021_b519.fasta")
original = read.fasta("original.txt")

mutations = data.frame(
  mutation = character(),
  pos_global = integer(),
  cambioCodon = character(),
  cambioAmino = character(),
  gen = character()
)

aminoacido = function(codon) {
  aminoacid = switch(
    codon, 
    "GCU" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "CGU" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGA" = "R", "AGG" = "R", "AAU" = "N", "AAC" = "N",
    "GAU" = "D", "GAC" = "D", "UGU" = "C", "UGC" = "C",
    "CAA" = "Q", "CAG" = "Q", "GAA" = "E", "GAG" = "E",
    "GGU" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G",
    "CAU" = "H", "CAC" = "H", "AUU" = "I", "AUC" = "I",
    "AUA" = "I", "UUA" = "L", "UUG" = "L", "CUU" = "L",
    "CUC" = "L", "CUA" = "L", "CUG" = "L", "AAA" = "K",
    "AAG" = "K", "AUG" = "M", "UUU" = "F", "UUC" = "F",
    "CCU" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "UCU" = "S", "UCC" = "S", "UCA" = "S", "UCG" = "S",
    "AGU" = "S", "AGC" = "S", "ACU" = "T", "ACC" = "T",
    "ACA" = "T", "ACG" = "T", "UGG" = "W", "UAU" = "Y",
    "UAC" = "Y", "GUU" = "V", "GUC" = "V", "GUA" = "V",
    "GUG" = "V"
  )
  return (aminoacid)
}


mutaciones = function(arnOriGenX, arnMexaGenX, inicioOri, gen){
  diff = which(arnOriGenX != arnMexaGenX) #obtener las posiciones donde hay nucleótidos diferentes
  for (x in diff){
    muta = paste(arnOriGenX[x],"to",arnMexaGenX[x], sep="") # formar un sting con el cambio de nucleótido
    inicioCodon = x - (x-1)%%3 # calcular el incio del codón
    posGlobal = inicioCodon + inicioOri #calcular la posición global
    numCodon = as.integer((x-1)/3+1) #calcular el número del codón
    indiceCodon = x - (x%%3) + 1
    codonOri = paste(arnOriGenX[indiceCodon], arnOriGenX[indiceCodon+1], arnOriGenX[indiceCodon+2], sep = "")
    codonMex = paste(arnMexaGenX[indiceCodon], arnMexaGenX[indiceCodon+1], arnMexaGenX[indiceCodon+2], sep = "")
    codon = paste(codonOri,"to",codonMex, sep="")
    aminoOri = aminoacido(codonOri)
    aminoMex = aminoacido(codonMex)
    amino = paste(aminoOri, numCodon, aminoMex, sep="")
    obs = list(muta, posGlobal, codon, amino, gen)
    mutations[nrow(mutations)+1, ] = obs
  }
  return(mutations)
}


for (k in seq(1, length(mx_2020))) {
  arnMexa2020 = as.vector(mx_2020[[k]])
  arnMexa2020[arnMexa2020=="t"] = "u"
  arnMexa2020 = toupper(arnMexa2020)
  for (j in seq(1,length(mx_2021))){
    if (k%%12 == j%%12){
      arnMexa2021 = as.vector(mx_2021[[j]])
      arnMexa2021[arnMexa2021 == "t"] = "u"
      arnMexa2021 = toupper(arnMexa2021)
      vect = j%%12
      if (vect == 0){
        vect = 12
      }
      if (j==2) next
      anotaciones = attr(original[[vect]], "Annot") 
      atributos = unlist(strsplit(anotaciones,"\\[|\\]|:|=|\\.|\\(")); 
      geneName = atributos[which(atributos=="gene")+1] 
      if (length(which(atributos=="join"))>0) inicioGen = as.integer(atributos[which(atributos=="join")+1]) 
      else inicioGen = as.integer(atributos[which(atributos=="location")+1]) 
      #cat ("------ gene:", geneName, "inicioGen:",inicioGen,"\n")
      mutations = mutaciones(arnMexa2020, arnMexa2021, inicioGen, geneName)
    }
  }
}


library(dplyr)
dfgraph = filter( #crea una especie de dataframe (tibble) con su propia estructura
  summarise(
    select(
      group_by(mutations, cambioAmino), #(df, columna), agrupa por columna de aminoacido
      mutation:gen
    ),
    mutation = first(mutation),
    cambioCodon = first(cambioCodon),
    gen = first(gen),
    Cuenta = n() #establece una cuenta para verificar el número de registros que tenemos
  ),
  Cuenta>100
)
mutations_S = subset(dfgraph, gen == "S")

#library(ggplot2)
#p = ggplot(mutations_S)
#p = p + aes(x=cambioAmino, y= Cuenta, fill=cambioAmino, label=..identity..)
#p = p + ggtitle("Frecuencia de mutaciones de sustitución en B.1.1.519 en el gen S")
#p = p + labs(x="Mutación", y="Frecuencia", fill="Mutación")
#p = p + geom_bar(stat = "identity") #stat = count es cuando ggplot hace el conteo
#p = p + geom_text(stat = "identity", vjust=0)
#p = p + theme_bw()

#