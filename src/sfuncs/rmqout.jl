# 
function rmqout()
        map(x->rm(x),glob("qoutjl"))
    end

    #ctlg(pwd(),"txt","weiss")

    """ 
    renamer - remove beginning char _   
    and replace it with C
    """
    