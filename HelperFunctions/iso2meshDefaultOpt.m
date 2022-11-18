function iso2meshOpt = iso2meshDefaultOpt()
    %% Default iso2mesh parameters
    iso2meshOpt.maxnode = 40000
    iso2meshOpt.radbound = 2
    iso2meshOpt.dofix=1
    iso2meshOpt.isovalues=0.5
    iso2meshOpt.method = 'cgalsurf' % 'cgalmesh' takes 3x more time (20 sec vs  60 sec)
    % otherwise, simply use: iso2meshOpt = 0.5;