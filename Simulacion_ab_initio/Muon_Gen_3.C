void disThet(){
    return;
}

void LandVav(){
    return ;
}

void Muon_Gen_3(){
    double half_sizex = 1.197 / 2; // cm
    double half_sizey = 1.587 / 2; // cm
    double half_sizez = 0.0725; // cm

    // Mapeos de las dimensiones 
    double stepxy = 0.0001;
    double stepz = 0.0001;

    int sizex =  (half_sizex * 2) /stepxy;
    int sizey =  (half_sizey * 2) /stepxy;
    int sizez =  (half_sizez * 2) /stepz;
    // std::cout << sizex  << std::endl;
    // std::cout << sizey  << std::endl;
    // std::cout << sizez  << std::endl;

    double mapx[sizex];
    double mapy[sizey];
    double mapz[sizez];


    // Mapeo de la CCD /// =======================> COMPLETADO
    double value =  - half_sizex;
    for (int index = 0; index <  sizex + 1 ; index++){
        if (index == 0){
            mapx[index] =  value;
            // std::cout <<  mapx[index]  << std::endl;
        }
        if (index > 0 ){
            value += stepxy;
            mapx[index] =  value;
            // std::cout <<  mapx[index]  << std::endl;
        }
    }

    value =  - half_sizey;
    for (int index = 0; index <  sizey + 2 ; index++){
        if (index == 0){
            mapy[index] =  value;
            // std::cout <<  mapy[index]  << std::endl;
        }
        if (index > 0 ){
            value += stepxy;
            mapy[index] =  value;
            // std::cout <<  mapy[index]  << std::endl;
        }
    }

    value =  - half_sizez;
    for (int index = 0; index <  sizez + 2 ; index++){
        if (index == 0){
            mapz[index] =  value;
            // std::cout <<  mapz[index]  << std::endl;
        }
        if (index > 0 ){
            value += stepz;
            mapz[index] =  value;
            // std::cout <<  mapz[index]  << std::endl;
        }
    }

    // double x = mapx[1];

    // std::cout << mapx  << std::endl;

    // if (half_sizez  mapz){

    // }

    return 0;
}