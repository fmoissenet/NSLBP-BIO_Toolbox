*Flow* est une extension de *pyCGM2* permettant de travailler avec les differents CGM
à partir d'un fichier de configuration optimisé pour l' analyse quantifiée de la marche.



Vous trouverez une documentation complete dans le repertoire *doc/build/html/index.html*


# Installation

Tout comme pyCGM2, l extension *flow* est un package **2.7**

La section installation de la documentation s'appuie sur l'environnement de programmation **Miniconda 3[64bits]**.
La documentation explique :

 * comment creer un environnement virtuel python2.7 - 32bits
 * comment installer *pyCGM2* et son extension *pyCGM2-Flow*


# Usage  

voici un usage concret de *flow* pour le CGM2_3

Considerons les données stockées dans *c:\\mes données\\Hannibal Lecter\\Session 1*. Toutes les c3d ont été labelisés et les evenements du cycle de marche identifiés.

Le traitement s'effectuera de la maniere suivante:

1. ouvrir anaconda prompt

2. activer votre environnement virtuel ( ici : pycgm2)

```
activate pycgm2
```

3. puis taper:

```
cd c:\\mes données\\Hannibal Lecter\\Session 1
pycgm2f-init.exe
pycgm2f-edit.exe -m CGM23
```

4. Editer le contenu de *CGM23.userSettings* et *emg.settings* puis taper :

```
pycgm2f-process.exe -m CGM23 -p
```
