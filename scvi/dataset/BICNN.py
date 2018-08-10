from scvi.dataset import LoomDataset


class MacoskoDataset(LoomDataset):
    def __init__(self, save_path='data/'):
        super(MacoskoDataset, self).__init__(filename='macosko.loom',
                                             save_path=save_path,
                                             url='https://github.com/YosefLab/scVI-data/raw/master/macosko.loom')


class RegevDataset(LoomDataset):
    def __init__(self, save_path='data/'):
        super(RegevDataset, self).__init__(filename='regev.loom',
                                           save_path=save_path,
                                           url='https://github.com/YosefLab/scVI-data/raw/master/regev.loom')


order = [['Pvalb low', 'Pvalb', 'Pvalb 1', 'Pvalb 2'],
         ['Pvalb Ex_1', 'Pvalb Ex_2', 'Pvalb Ex'],
         ['Pvalb Astro_1', 'Pvalb Astro_2'],
         ['L2/3 IT Astro', 'L2/3 IT Macc1', 'L2/3 IT Sla_Astro', 'L2/3 IT', 'L2/3 IT Sla', 'L2/3 IT Sla_Inh'],
         ['Sst Tac2', 'Sst Myh8', 'Sst Etv1', 'Sst Chodl', 'Sst'],
         ['L5 PT_2', 'L5 PT IT', 'L5 PT_1'],
         ['L5 IT Tcap_1_3', 'L5 IT Tcap_2', 'L5 IT Tcap_Astro', 'L5 IT Tcap_1', 'L5 IT Tcap_L2/3',
          'L5 IT Tcap_Foxp2', 'L5 IT Tcap_3'],
         ['L5 IT Aldh1a7_2', 'L5 IT Aldh1a7', 'L5 IT Aldh1a7_1'],
         ['L5 NP', 'L5 NP Slc17a8'],
         ['L6 IT Car3', 'L6 CT Olig', 'L6 IT Maf', 'L6 IT Ntn5 Mgp', 'L6 IT Ntn5 Inpp4b'],
         ['L6 CT Nxph2', 'L6 CT Astro', 'L6 CT', 'L6 CT Grp'],
         ['L6b', 'L6b F2r'],
         ['Lamp5 Sncg', 'Lamp5 Egln3', 'Lamp5 Slc35d3'],
         ['Vip Rspo4', 'Vip Serpinf1', 'Vip'],
         ['Astro Ex', 'Astro Aqp4'],
         ['OPC Pdgfra'],
         ['VLMC Osr1'],
         ['Oligo Enpp6_1', 'Oligo Enpp6_2', 'Oligo Opalin'],
         ['Sncg Ptprk'],
         ['Endo Slc38a5', 'Endo Slc38a5_Peri_2', 'Endo Slc38a5_Peri_1']]

order = [[key for key in key_color_group] for key_color_group in order]
ravel = [key for key_color in order for key in key_color]
labels_group_array = [len(key) * [group] for group, key in enumerate(order)]
labels_groups = []
for group in labels_group_array:
    labels_groups += group
