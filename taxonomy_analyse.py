import pandas as pd

df = pd.read_excel("/mnt/c/Users/jpsfr/Downloads/Taxonomy_venn.xlsx", sheet_name="venn")

metaproteomes = df["Metaproteomes"].tolist()
metaproteomes = [metaproteome.strip() for metaproteome in metaproteomes]
metaproteomes = sorted(list(set(metaproteomes)))

metagenomes = df["Metagenomes"].tolist()
metagenomes = [metagenome.strip() for metagenome in metagenomes if not isinstance(metagenome, float)]
metagenomes = list(set(metagenomes))


common_elements = list(filter(lambda x: x in metaproteomes, metagenomes))