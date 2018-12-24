double getCost_mybin(int i){
    const char *command01 = catStrIntStr(tempFileAddress, "temp_", i, ".pdb");
    const char *command02 = catStrIntStr(mybinAddress, "temp_", i, ".pdb");
    const char *command0 = catStrStr("cp ", command01, " ", command02);
    system(command0);

    double ans=0;
    while (ans==0){
        int lineNum=0;
        char buffer[10000];
        //system("/home/ws/GL/mybin/calcquarkenergy /home/ws/GL/mybin/ /home/ws/GL/mybin/temp1.pdb>>QUACKout.txt");
        const char *str = catStrStr(mybinAddress, "calcquarkenergy ", mybinAddress, " ",
                                    mybinAddress, "temp_");
        const char *command1 = catStrIntStr(str, i, ".pdb>>");
        const char *command2 = catStrIntStr(QUACKoutFileAddress, "QUACKout", i, ".txt");
        const char *command = catStrStr(command1, command2);
        system(command);
        ifstream infile;
        openFile(command2, infile);
        while (infile.getline(buffer, 10000)) lineNum++;
        infile.close();
        openFile(command2, infile);
        for (int i=0; i<lineNum-1; i++) infile.getline(buffer, 10000);
        infile >> ans;
        infile.close();
        removeFile(command2);
        delete command;  delete command1; delete command2;
    }
    removeFile(command02);

    delete command0;  delete command01; delete command02;

    return ans;
}
