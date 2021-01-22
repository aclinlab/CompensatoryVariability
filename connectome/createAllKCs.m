for i=0:1926
    kc = KCskel.createKC(i);
    kc = kc.getDistBetweenPNInputs;
    kc = kc.getDistPostPedToPNInputs;
    kc.saveKCskel;
end