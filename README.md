# ProjetoFinal-WBDSLA
Esse repositório contém o código e desenvolvimento do projeto final do WBDS LA Camp, um curso de treinamento intensivo gratuito para estudantes de graduação e pós-graduação para aprender bioinformática e ciência de dados.

## Resumo

Construí um modelo de Machine Learning com o objetivo de classificar amostras de Astrocitoma e Glioblastoma baseado em seus dados de metilação de DNA (beta-values). O modelo escolhido foi o Random Forest e o projeto foi desenvolvido integralmente na linguaguem de programação Python. <br><br>
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

Com isso, foi possível observar o número total de amostras (N=145). <br><br>
A partir desta análise exploratória, foi feita a deleção de 20 amostras que não continham dados de metilação de DNA, diminuindo o espaço amostral para 120. <br><br>




