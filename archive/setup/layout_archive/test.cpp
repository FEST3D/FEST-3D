Grid[ ]; // array of grids
i = 0;
while(Grid[i])
{
    j = 0;
    while(Grid[j])
    {
        if(Grid[i].IMIN_FACE.coordinates == Grid[j].IMAX_FACE.coordinates)
        {
            Populate_interface_info(Grid[i],Grid[j],IMIN);
        }
        if(Grid[i].JMIN_FACE.coordinates == Grid[j].JMAX_FACE.coordinates)
        {
            Populate_interface_info(Grid[i],Grid[j],JMIN);
        }
        if(Grid[i].KMIN_FACE.coordinates == Grid[j].KMAX_FACE.coordinates)
        {
            Populate_interface_info(Grid[i],Grid[j],KMIN);
        }
        j++;
    }
    i++;
}
