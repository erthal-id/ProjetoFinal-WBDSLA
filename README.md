# ProjetoFinal-WBDSLA
Esse repositório contém o código e desenvolvimento do projeto final do WBDS LA Camp, um curso de treinamento intensivo gratuito para estudantes de graduação e pós-graduação para aprender bioinformática e ciência de dados.

## Resumo

Construí um modelo de Machine Learning com o objetivo de classificar amostras de Astrocitoma e Glioblastoma baseado em seus dados de metilação de DNA (beta-values) 27k. O modelo escolhido foi o Random Forest e o projeto foi desenvolvido integralmente na linguaguem de programação Python. <br><br>
Realizei a aquisição, pré-processamento e tratamento dos dados, além do desenvolvimento do modelo. <br><br>
Neste repositório, há o Jupyter Notebook e o arquivo .py, ambos com o código para que seja possível a reprodução do projeto.

## Instalação

No bloco de código abaixo estão os pacotes, classes e funções utilizados ao longo do projeto com uma breve descrição de suas finalidades

```python 
import GEOparse #Download dos dados do NCBI/GEO
import pandas as pd #Tratamento dos dados
import seaborn as sns #Gráficos
from matplotlib import pyplot as plt #Gráficos
from sklearn.metrics import accuracy_score #Cálculo da acurácia
from sklearn.ensemble import RandomForestClassifier #Modelo utilizado
from sklearn.model_selection import train_test_split #Partição do dataset
from sklearn.feature_selection import SelectKBest, chi2 #Seleção de features
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay #Construção da matriz confusão
```

## Execução

Meu projeto não há arquivos de entrada e saída que necessariamente precisam estar na máquina de quem está o executando. <br><br>
O dado deste projeto pode ser encontrado no banco de dados GEO (Gene Expression Omnibus), do NCBI (National Center for Biotechnology Information), através do ID: GSE53229. <br><br>
O ID foi utilizado para a aquisição desse dado para o ambiente de execução, como foi feito a seguir: 

```python 
gse = GEOparse.get_GEO(geo="GSE53229", destdir="./")
```

E a partir do objeto ``` gse ``` é possível obter um dataframe com diversas informações acerca do dataset, a partir do comando:
```python 
phenodata = gse.phenotype_data 
```

Com isso, foi possível observar o número total de amostras (N=145).

### Seleção do espaço amostral

A partir desta análise exploratória, foi feita a deleção de 20 amostras que não continham dados de metilação de DNA, diminuindo o espaço amostral para 120. <br><br>
Dentre essas 120 amostras, foram filtradas as amostras da classe "Human Astrocytoma" e "Human Glioblastoma", que continham 43 e 37 amostras, respectivamente. <br><br>
Logo, o espaço amostral desse projeto foi de 80 amostras, divididas entre amostras de Astrocitoma e Glioblastoma. <br><br>

```python 
DNAmeth_phenodata = phenodata[phenodata['type']=='genomic']
astrocitoma = list(DNAmeth_phenodata.index[DNAmeth_phenodata["source_name_ch1"] == "Human astrocytoma"])
glioblastoma = list(DNAmeth_phenodata.index[DNAmeth_phenodata["source_name_ch1"] == "Human glioblastoma"])
samples = astrocitoma + glioblastoma
```

O objeto ``` samples ``` é uma lista com o ID das 80 amostras que foram utilizadas. Essa lista será útil para o próximo passo.

### Contrução da matriz de níveis de metilação de DNA das amostras

Foi construída a função ``` formar_dataframe() ```, em que é acessado o nível de metilação em 27.578 locais/probes no DNA para cada amostra, e que tem como argumento apenas uma lista com os IDs das amostras:

```python 
def formar_dataframe(lista_de_amostras):

  df = pd.DataFrame(gse.gsms[lista_de_amostras[0]].table[["ID_REF","VALUE"]]) #acessa o nível de metilação da amostra em questão
  df.rename(columns={"VALUE":lista_de_amostras[0]}, inplace = True) #renomeia coluna "VALUE" para o nome da amostra
  df = df.set_index("ID_REF") #transformo a coluna "ID_REF" em rownames
  df = df.T #faz a transposta do dataframe

  for i in range(1,len(lista_de_amostras)):
    df2 = pd.DataFrame(gse.gsms[lista_de_amostras[i]].table[["ID_REF","VALUE"]])
    df2.rename(columns={"VALUE":lista_de_amostras[i]}, inplace = True)
    df2 = df2.set_index("ID_REF")
    df2 = df2.T 
    if list(df.columns) == list(df2.columns): #confere se a ordem das probes é a mesma 
      df = pd.concat([df, df2]) #concatena os dataframes iterativamente
    else:
      df2 = df2[df.columns] #ordena caso não seja
      df = pd.concat([df, df2])  
  
  return df
  ```
  
  E com isso, foi possível construir o dataframe principal deste projeto, que é a matriz de metilação de DNA, com dimensão de 80 X 27.579, sendo uma coluna adicionada posteriormente para rotular os tipos das amostras entre Astrocitoma e Glioblastoma:
  
  ```python
  
  methylation_matrix = formar_dataframe(samples)
  label = list((["1"] * 43) + ["2"] * 37) #1: Astrocitoma, 2: Glioblastoma
  methylation_matrix["label"] = label
  ```
  
  O dataframe gerado fica neste seguinte formato:
  
| ID_REF | cg07259382 | cg17271365 | ... | label |
|------------|----------|----------|-----|---|
| GSM1287993 | 0.379213 | 0.557644 | ... | 1 |
| GSM1287994 | 0.514338 | 0.445652 | ... | 1 |
| GSM1288116 | 0.807990 | 0.429600 | ... | 2 |
| GSM1288117 | 0.220680 | 0.250000 | ... | 2 |

80 rows X 27579 columns




