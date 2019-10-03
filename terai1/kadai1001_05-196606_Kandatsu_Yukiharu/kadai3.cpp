/*kadai3 05-195506 Yukiharu Kandatsu*/
#include<iostream>
#include<fstream>
#include<string>

using namespace std;

#define inf 1000000

//function for compare two letter.
int s(char x, char y, int match, int mismatch){
    if(x == y) return match;
    else return mismatch;
}

void show(int sizex, int sizey, int*table){
    ofstream ost("xxx.txt");
    for(int j = 0; j < sizey + 1; j++){
        for(int i = 0; i < sizex + 1; i++){
            ost << table[j*(sizex+1)+i] << " ";
        }
        ost << endl;
    }
    ost.close();
}

void show2(int sizex, int sizey, int*table){
    for(int j = 0; j < sizey + 1; j++){
        for(int i = 0; i < sizex + 1; i++){
            cout << table[j*(sizex+1)+i] << " ";
        }
        cout << endl;
    }
}

int main(void){
    //begin file input.
    ifstream ist("alignment.txt");
    string seqX; string seqY;

    string checker;
    getline(ist, checker); //delete first line.
    while(!ist.eof()){
        getline(ist, checker);
        if(checker[0] == '>'){
            break;
        }
        else seqX += checker;
    }
    while(!ist.eof()){
        getline(ist, checker);
        seqY += checker;
    }
    //end file input.

    //begin: define variables and DP table and TraceBack table.
    int sizex = seqX.length(); int sizey = seqY.length();

    int match = 1; int mismatch = -1; int d = 5; int e = 1;

    int* dptM = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* dptX = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* dptY = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* bdptM = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* bdptX = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* bdptY = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    /*
    int* tbtM = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* tbtX = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    int* tbtY = (int*)malloc(sizeof(int) * (sizex + 1) * (sizey + 1));
    */
    int* A = (int*)malloc(sizeof(int) * (sizex) * (sizey));
    //end: define variables and DP table and TraceBack table.

    //begin: initialize
    for(int i = 1; i < sizex + 1; i++){
        dptM[i] = -inf;
        dptX[i] = -d - (i - 1)*e;
        dptY[i] = -inf;
        //tbtX[i] = 1;
    }
    for(int j = 1; j < sizey + 1; j++){
        dptM[j * (sizex + 1)] = -inf;
        dptX[j * (sizex + 1)] = -inf;
        dptY[j * (sizex + 1)] = -d - (j - 1)*e;
        //tbtY[j*(sizex+1)] = 2;
    }
    dptM[0] = 0; dptX[0] = -inf; dptY[0] = -inf;

    //tbtM[0] = -1; tbtX[0] = -1; tbtY[0] = -1;
    //end : initialize 

    //begin: fill table
    for(int i = 1; i < sizex + 1; i++){
        for(int j = 1; j < sizey + 1; j++){
            //begin: M(i, j)
            int candidates[3] = {-inf};
            candidates[0] = dptM[(j-1)*(sizex+1)+i-1] + s(seqX[i-1], seqY[j-1], match, mismatch);
            candidates[1] = dptX[(j-1)*(sizex+1)+i-1] + s(seqX[i-1], seqY[j-1], match, mismatch);
            candidates[2] = dptY[(j-1)*(sizex+1)+i-1] + s(seqX[i-1], seqY[j-1], match, mismatch);
            int biggest = candidates[0];
            int ori = 0;
            for(int i = 1; i < 3; i++){
                if(candidates[i] > biggest){
                    biggest = candidates[i];
                    ori = i;
                }
            }
            dptM[j*(sizex+1)+i] = biggest;
            //tbtM[j*(sizex+1)+i] = ori;
            //end: M(i,j)

            //begin: X(i, j)
            candidates[0] = -inf; candidates[1] = -inf;
            candidates[0] = dptM[j*(sizex+1)+i-1] + - d;
            candidates[1] = dptX[j*(sizex+1)+i-1] + - e;
            ori = 0;
            if(candidates[1] > candidates[0]) ori = 1;
            dptX[j*(sizex+1)+i] = candidates[ori];
            //tbtX[j*(sizex+1)+i] = ori;
            //end: X(i, j)

            //begin: Y(i, j)
            candidates[0] = -inf; candidates[2] = -inf;
            candidates[0] = dptM[(j-1)*(sizex+1)+i] + - d;
            candidates[2] = dptY[(j-1)*(sizex+1)+i] + - e;
            ori = 0;
            if(candidates[2] > candidates[0]) ori = 2;
            dptY[j*(sizex+1)+i] = candidates[ori];
            //tbtY[j*(sizex+1)+i] = ori;
            //end: Y(i, j)
        }
    }
    //end: fill table

    //begin: initialize for back
    for(int i = sizex - 1; i > -1; i--){
        bdptM[(sizey)*(sizex + 1) + i] = -d - (sizex - i - 1)*e;
        bdptX[(sizey)*(sizex + 1) + i] = -(sizex - i) * e;
        bdptY[(sizey)*(sizex + 1) + i] = -inf;
    }
    for(int j = sizey - 1; j > -1; j--){
        bdptM[j * (sizex + 1) + sizex] = -d - (sizey - j - 1)*e;
        bdptX[j * (sizex + 1) + sizex] = -inf;
        bdptY[j * (sizex + 1) + sizex] = - (sizey - j)*e;
    }
    bdptM[(sizey)*(sizex + 1) + (sizex)] = 0; bdptX[(sizey)*(sizex + 1) + (sizex)] = 0; bdptY[(sizey)*(sizex + 1) + (sizex)] = 0;

    //tbtM[0] = -1; tbtX[0] = -1; tbtY[0] = -1;
    //end : initialize for back

    //begin: fill table for back
    for(int i = sizex - 1; i > -1; i--){
        for(int j = sizey - 1; j > -1; j--){
            //begin: M(i, j)
            int candidates[3] = {-inf};
            candidates[0] = bdptM[(j+1)*(sizex+1)+i+1] + s(seqX[i], seqY[j], match, mismatch);
            candidates[1] = bdptX[(j)*(sizex+1)+i+1] - d;
            candidates[2] = bdptY[(j)*(sizex+1)+i+1] - d;
            int biggest = candidates[0];
            int ori = 0;
            for(int i = 1; i < 3; i++){
                if(candidates[i] > biggest){
                    biggest = candidates[i];
                    ori = i;
                }
            }
            bdptM[j*(sizex+1)+i] = biggest;
            //tbtM[j*(sizex+1)+i] = ori;
            //end: M(i,j)

            //begin: X(i, j)
            candidates[0] = -inf; candidates[1] = -inf;
            candidates[0] = bdptM[(j + 1)*(sizex+1)+i+1] + s(seqX[i], seqY[j], match, mismatch);
            candidates[1] = bdptX[j*(sizex+1)+i+1] - e;
            ori = 0;
            if(candidates[1] > candidates[0]) ori = 1;
            bdptX[j*(sizex+1)+i] = candidates[ori];
            //tbtX[j*(sizex+1)+i] = ori;
            //end: X(i, j)

            //begin: Y(i, j)
            candidates[0] = -inf; candidates[2] = -inf;
            candidates[0] = bdptM[(j+1)*(sizex+1)+i+1] + s(seqX[i], seqY[j], match, mismatch);
            candidates[2] = bdptY[(j+1)*(sizex+1)+i] - e;
            ori = 0;
            if(candidates[2] > candidates[0]) ori = 2;
            bdptY[j*(sizex+1)+i] = candidates[ori];
            //tbtY[j*(sizex+1)+i] = ori;
            //end: Y(i, j)
        }
    }
    //end: fill table for back

    //begin: fill answer
    for(int i = 0; i < sizex; i++){
        for(int j = 0; j < sizey; j++){
            A[(j)*(sizex)+i] = dptM[(j+1)*(sizex + 1)+i + 1] + bdptM[(j + 1)*(sizex+1)+i + 1];
        }
    }
    //end: fill answer

    /* /begin traceback
    string ansX; string ansY; string hy ="-";
    int final[3] = {0};
    final[0] = dptM[(sizey+1)*(sizex+1)-1];
    final[1] = dptX[(sizey+1)*(sizex+1)-1];
    final[2] = dptY[(sizey+1)*(sizex+1)-1];
    int ori = 0;
    for(int i = 1; i < 3; i++){
        if(final[i] > final[ori]) ori = i;
    }

    int curr = ori;
    int i = sizex - 1; int j = sizey - 1;
    int count = 0;
    while(curr != -1 && (i >= 0 || j >= 0)){
        if(curr == 0){
            ansX.insert(0, seqX, i, 1); ansY.insert(0, seqY, j, 1);
            curr = tbtM[(j+1)*(sizex+1)+(i+1)];
            i--; j--;
        }
        else if(curr == 1){
            ansX.insert(0, seqX, i, 1); ansY.insert(0, hy);
            curr = tbtX[(j+1)*(sizex+1)+(i+1)];
            i--;
        }
        else if(curr == 2){
            ansX.insert(0, hy); ansY.insert(0, seqY, j, 1);
            curr = tbtY[(j+1)*(sizex+1)+(i+1)];
            j--;
        }
        count++;
    }

    int p = 40; //num of letter to show in one line
    for(int k = 0; k*p < count; k++){
        cout << "X : ";
        for(int l = 0; l < p; l++){
            if(ansX[k*p + l] == 0) break;
            cout << ansX[k*p + l];
        }
        cout << endl;
        cout << "Y : ";
        for(int l = 0; l < p; l++){
            if(ansY[k*p + l] == 0) break;
            cout << ansY[k*p + l];
        }
        cout << endl << endl;
    }
    //end traceback
    */

    /*
    
    show(sizex, sizey, dptX);
    cout << endl;
    show(sizex, sizey, dptY);
    cout << endl << endl;
    show(sizex, sizey, bdptM);
    cout << endl;
    show(sizex, sizey, bdptX);
    cout << endl;
    show(sizex, sizey, bdptY);
    cout << endl << endl;
    show(sizex, sizey, bdptM);
    cout << endl;
    
    show(sizex, sizey, dptM);
    cout << endl;
    show2(sizex, sizey, dptM);
    cout << endl;
    show2(sizex, sizey, bdptM);
    cout << endl;
    
    
    show(sizex, sizey, tbtM);
    cout << endl;
    show(sizex, sizey, tbtX);
    cout << endl;
    show(sizex, sizey, tbtY);
    */
    
    show(sizex - 1, sizey - 1, A);
    cout << endl;
    return 0;
}