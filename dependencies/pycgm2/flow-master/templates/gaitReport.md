% Compte rendu d'Analyse de la marche
% **IPP :** {{ data["infos"]["Ipp"] }}
% **Numero Visite :** {{ data["infos"]["SessionNumber"] }}

**date :** {{ data["infos"]["Date"] }}
**Age :** {{ data["infos"]["Age"] }} ans


## Analyses

**modèle biomécanique :** {{ data["infos"]["Model"] }}


{% for anaIt in  data['Report'].Analyses %}
### {{ anaIt["section"] }}
    {% for it in  anaIt["pages"] %}
#### {{ anaIt["pageTitles"][loop.index-1]}}
![]({{it}}){width=14cm}
    {% endfor %}
{% endfor %}

## Comparaisons

{% for compIt in  data['Report'].Comparisons %}
### {{ compIt["section"] }}
    {% for it in  compIt["pages"] %}
#### {{ compIt["pageTitles"][loop.index-1]}}
![]({{it}}){width=14cm}
    {% endfor %}
{% endfor %}
