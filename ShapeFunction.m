
function [N, DNL]=ShapeFunction(Ng,elementType)


switch elementType

    
    case '1D'

N=[0.5*(1-Ng),0.5*(1+Ng)];
DNL=[-0.5 0.5];

    
    case 'TRI3'

ksi=Ng(1);eta=Ng(2);        
N=[1-ksi-eta,ksi,eta];
DNL_ksi=[-1, 1,0];
DNL_eta=[-1,0,1];
DNL=[DNL_ksi;DNL_eta];

    case 'QUAD4'

ksi=Ng(1);eta=Ng(2);
N=0.25*[(1-ksi)*(1-eta),(1+ksi)*(1-eta),(1+ksi)*(1+eta),(1-ksi)*(1+eta)];

DNL_ksi=0.25*[-(1-eta) (1-eta) (1+eta) -(1+eta)];
DNL_eta=0.25*[-(1-ksi) -(1+ksi) (1+ksi) (1-ksi)];
DNL=[DNL_ksi;DNL_eta];

    otherwise
        error('Unknown element type')
end

end