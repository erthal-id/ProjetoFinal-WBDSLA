# ProjetoFinal-WBDSLA
Esse repositório contém o código e desenvolvimento do projeto final do WBDS LA Camp, um curso de treinamento intensivo gratuito para estudantes de graduação e pós-graduação para aprender bioinformática e ciência de dados.

## Resumo

Construí um modelo de Machine Learning com o objetivo de classificar amostras de Astrocitoma e Glioblastoma baseado em seus dados de metilação de DNA (beta-values) 27k. O modelo escolhido foi o Random Forest e o projeto foi desenvolvido integralmente na linguaguem de programação Python. <br><br>
Realizei a aquisição, pré-processamento e tratamento dos dados, além do desenvolvimento do modelo. <br><br>
Neste repositório, há o Jupyter Notebook e o arquivo .py, ambos com o código para que seja possível a reprodução do projeto.

## Justificativa e Objetivo
A Organização Mundial da Saúde já incentiva as análises de eventos moleculares, como a  metilação de DNA, para a classificação de tumores do sistema nervoso central. <br><br>
Meu objetivo com este projeto final foi aplicar o que foi visto ao longo do Camp - Análise de dados e Ciência de Dados, à análise de metilação de DNA para tentar classificar amostras de Astrocitoma e Glioblastoma, construindo um algoritmo supervisionado Random Forest.

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
| ... | ... | ... | ... | ... |
| GSM1288116 | 0.220680 | 0.250000 | ... | 2 |
| GSM1288117 | 0.220680 | 0.250000 | ... | 2 |

80 rows X 27579 columns <br><br>

Como etapa de tratamento dos dados, foram removidos os dados faltantes, pois eram poucos: 

```python 
methylation_matrix = methylation_matrix.dropna(axis = "columns") #Remove colunas com NAs
methylation_matrix.shape #Nova dimensão da matriz de metilação de DNA: 27.274 probes (col) X 80 samples (row)
```
 ### Seleção de features
 
Como há um número muito alto do que seriam as features (27.274), foi feita uma seleção de features baseada em filtro para diminuir esse número antes da construção do modelo. Essa etapa foi realizada utilizando o teste qui-quadrado, que seleciona as 1000 melhores features:

```python 
#Separando as features da label
X = methylation_matrix.iloc[:, :-1]
y = methylation_matrix.iloc[:, -1]

#Seleção de features baseada em filtro
selector = SelectKBest(chi2, k=1000) #Seleciona as 1000 melhores features usando o teste qui-quadrado
X_new = selector.fit_transform(X, y) #seleciona as colunas no dataframe original

#Obtendo as colunas selecionadas
mask = selector.get_support() #boolean array com as colunas selecionadas
new_features = X.columns[mask] #array com os nomes das colunas selecionadas

#Criando o novo dataframe apenas com as colunas selecionadas e o rótulo
filtered_methylation_matrix = methylation_matrix[list(new_features)]
filtered_methylation_matrix["label"] = label
```

### Construção do modelo 

Foi feita a partição do dataset em 75% treino e 25% teste e, dessa forma, foi construído o modelo com os parâmetros padrões:

```python
SEED = 42

X = filtered_methylation_matrix.drop(columns=["label"])
y = filtered_methylation_matrix.label

#25% para teste
Xtrain, Xtest, ytrain, ytest = train_test_split(X, y, test_size = 0.25, stratify = y, random_state = SEED)
#Treinamento: 60, Teste: 20

#Treino e teste com os parâmetros padrões
model = RandomForestClassifier()
model.fit(Xtrain,ytrain)
acc = accuracy_score(ytest, model.predict(Xtest))
```
A acurácia obtida, armazenada em ``` acc ``` foi de 0.9 (90%).<br><br>
Como podemos observar na seguinte matriz confusão, houveram 11 Verdadeiros Positivos, 0 Falsos Negativos, 2 Falsos Positivos e 7 Verdadeiros Negativos: <br>
![image](https://user-images.githubusercontent.com/49324017/226037713-64ae229c-6365-4292-8f8f-4e20bff2adb7.png)

<br>
Código para a matriz de confusão:

```python 
ConfusionMatrixDisplay(confusion_matrix=confusion_matrix(ytest, model.predict(Xtest)), display_labels=model.classes_).plot()
plt.title("Matriz confusão\n ")
plt.show()
```

## Conclusões

Como é possível observar, o modelo obteve uma taxa alta de acurácia. Isso pode indicar, sim, que foi um bom classificador, no entanto eu me aprofundaria nas questões de Etapa de Validaçao e Ajuste de Hiperparâmetros para me certificar que não houve sobreajuste do modelo no meu conjunto de treino. Também não tenho certeza acerca da forma que escolhi realizar a seleção de features, visto que eu não tenho muito conhecimento sobre isso hoje. <br><br>
Dessa forma, posso dizer que o projeto foi muito importante para a fixação do conteúdo e treinamento de tudo que aprendi e, por mais que não verifiquei as questões citadas acima, fico satisfeita com o aprendizado.





