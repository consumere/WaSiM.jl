global tot = 0;
wd=pwd();
for (root, dirs, files) in walkdir(wd)
    s = [joinpath(root, file) for file in files]
    fzs = [(file, filesize(file) / 2^20) for file in s]
    [println("summing up $d ...") for d in dirs]
    global tot += sum(map(x->x[2], fzs))
end
tot=round(tot; digits=3);
printstyled("Size of $wd and its subdirectories: $tot MB\n", color=:green)