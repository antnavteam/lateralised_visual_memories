function poloarplot2rose (thref,activation)

data_rose = []; 

k=0;
    for i=thref %loop through all angles
    k=k+1;

    nb_data = round(activation(k)*1000);
    a=ones(1,nb_data).* i;

    data_rose= [data_rose a];
    end

rose(data_rose,thref)

