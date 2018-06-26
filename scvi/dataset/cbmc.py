from .csv import CsvDataset


class CbmcDataset(CsvDataset):
    def __init__(self, save_path='data/'):
        super(CbmcDataset, self).__init__(filename='GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz',
                                          save_path=save_path,
                                          compression='gzip',
                                          url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format="
                                              "file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz")
