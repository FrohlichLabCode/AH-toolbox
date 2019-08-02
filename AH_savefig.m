function AH_savefig(fig, directory, saveName)

AH_mkdir(directory)
savefig(fig, [directory saveName '.fig'],'compact');
saveas(fig, [directory saveName '.png']);