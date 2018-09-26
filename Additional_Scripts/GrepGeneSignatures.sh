#"cd4_t_helper", "regulatory_t", "naive_t", "memory_t", "cytotoxic_t", "naive_cytotoxic","b_cells", "cd34", "cd56_nk", "cd14_monocytes"

#cd4_t_helper - naive_cytotoxic
grep -G ".*CD8.*UP\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*CD4.*UP\thttp.*"| grep -v ".*DC.*UP\thttp.*"| grep -v ".*KO.*UP\thttp.*" |  grep -v ".*PATIENT.*UP\thttp.*"|  grep -v ".*COSTIM.*UP\thttp.*"| grep -v ".*INT.*UP\thttp.*" > CD4CD8.UP.gmt
grep -G ".*CD8.*DN\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*CD4.*DN\thttp.*"| grep -v ".*DC.*DN\thttp.*"| grep -v ".*KO.*DN\thttp.*" |  grep -v ".*PATIENT.*DN\thttp.*"|  grep -v ".*COSTIM.*DN\thttp.*"| grep -v ".*INT.*DN\thttp.*" > CD4CD8.DN.gmt


# naive_cytotoxic - b_cells 
grep -G ".*CD8.*DN\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*BCELL.*DN\thttp.*" > CD8BCELL.DN.gmt
grep -G ".*CD8.*UP\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*BCELL.*UP\thttp.*" > CD8BCELL.UP.gmt

# cd4_t_helper - b_cells 
grep -G ".*CD4.*DN\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*BCELL.*DN\thttp.*" | grep -v ".*STIM.*DN\thttp.*" |grep -v ".*LUPUS.*DN\thttp.*"|grep -v ".*PRO.*DN\thttp.*" |grep -v ".*MEMORY.*DN\thttp.*"  > CD4BCELL.DN.gmt
grep -G ".*CD4.*UP\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*BCELL.*UP\thttp.*" | grep -v ".*STIM.*UP\thttp.*" |grep -v ".*LUPUS.*UP\thttp.*"|grep -v ".*PRO.*UP\thttp.*" |grep -v ".*MEMORY.*UP\thttp.*"  > CD4BCELL.UP.gmt

# naive_cytotoxic - NK 
grep -G ".*CD8.*DN\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*NK.*DN\thttp.*" |grep -v ".*NKT.*DN\thttp.*" |grep -v ".*DC.*DN\thttp.*"|grep -v ".*INFECTION.*DN\thttp.*" > CD8NK.DN.gmt
grep -G ".*CD8.*UP\thttp.*" c7.all.v6.2.symbols.gmt|grep ".*NK.*UP\thttp.*" |grep -v ".*NKT.*UP\thttp.*" |grep -v ".*DC.*UP\thttp.*"|grep -v ".*INFECTION.*UP\thttp.*" > CD8NK.UP.gmt

