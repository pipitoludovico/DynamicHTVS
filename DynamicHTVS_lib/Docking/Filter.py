class Filter:
    def __init__(self, df, idCol):
        self.idCol = idCol
        self.df = df.drop('smiles', axis=1)
        self.dfclean = self.df.drop_duplicates(subset=['CanonicalSmiles'])

    def getFiltered(self):
        return self.dfclean.loc[:, (self.idCol, 'CanonicalSmiles')]
