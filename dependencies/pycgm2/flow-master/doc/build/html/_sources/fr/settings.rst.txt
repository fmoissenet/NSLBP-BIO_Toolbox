.. _userSettings:

User Settings
===============

.. contents::
    :depth: 3


Chaque CGMi possede son fichier **userSettings**. Le but de ce fichier est double. Il vise à :

 - renseigner les données d entrées du modele ( ex : les donnees anthropometriques, les options de calibration)
 - renseigner les acquisitions par condition d'analyse.

Ce fichier de configuration se veut le plus concis possible de maniere a rapidement dsitinguer le traitement ayant été réalisé.

Vous trouverez 2 niveaux de settings ( basic et Guru).

**Je recommande de maitriser les fichiers Basic**. Les fichiers Guru necessitent une connaissance profonde des CGM2.i

Basic userSettings
------------------------

Voici exemple de fichier *basic* settings. Il est ecrit en yaml.

.. literalinclude:: C:/Users/fleboeuf/Documents/Programmation/pyCGM2/pyCGM2-extensions/pyCGM2_flow/doc/details_basicSettings.yaml
    :language: yaml

Il est important de noter les relations "1-n" entre :

- les elements de *Calibration* et les differents *Trials* de la section *Fitting*. *UNE* calibration peut etre appliquée a *N* trials de marche

-  les elements de *Condition* et les differents *Trials* des la section *Fitting* et *emg*. UNE condition de marche pointe sur *N* trials de marche



Guru userSettings
------------------------

TODO
