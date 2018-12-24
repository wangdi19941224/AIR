//
// Created by advanced on 17-1-16.
//

#include "MOPSO.h"

void sortRepByAngle(myRep &rep, ANGLE *angle){
    int len = static_cast<int>(rep.size());
    myRep::iterator it ;
    POINT point[len+100];
    int k=0;
    for (it = rep.begin(); it != rep.end(); it++){
        point[k].x = it->Cost[0];
        point[k].y = it->Cost[1];
        point[k].id = k;
        k++;
    }
    pointQsort(0, len-1, point);
    for (int i=0; i<len; i++){
        angle[i].id = point[i].id;
    }
    for (int i=1; i<len-1; i++){
        getAngle(i, point, angle);
    }
    angle[0].value     = -PI;
    angle[len-1].value = -2*PI;
    angleQsort(0, len - 1, angle);

    for (int i = 0; i < len; i++){
        sortAns[i] = angle[i].id;
    }
}

void getAngle(int i, const POINT *point, ANGLE *angle){
    VECTOR AB(point[i-1].x-point[i].x, point[i-1].y-point[i].y);
    VECTOR AC(point[i+1].x-point[i].x, point[i+1].y-point[i].y);

    double DOT = vectorDot(AB, AC);
    double l1 = length(AB), l2 = length(AC);
    double cosTheta = DOT / (l1 * l2);
    angle[i].value = acos(cosTheta);
    if (pdPlusPi(point[i-1], point[i], point[i+1]))
        angle[i].value = 2*PI - angle[i].value;
}

bool pdPlusPi(POINT first, POINT second, POINT third) {
    VECTOR BA(first.x - second.x, first.y - second.y), BC(third.x - second.x, third.y - second.y);
    if (cross(BA,BC) < 0) return 1;
    else return 0;
}

double vectorDot(VECTOR u, VECTOR v){
    return u.x * v.x + u.y * v.y;
}

double cross(VECTOR u, VECTOR v){
    return u.x * v.y - v.x * u.y ;              //  x1*y2 - x2*y1
}

double length(VECTOR u){
    return sqrt(u.x*u.x + u.y*u.y);
}

void angleQsort(int l, int r, ANGLE *a){
    if (l>=r)
        return;

    int i=l, j=r;
    double holeValue=a[i].value;
    int holeId = a[i].id;
    while (i<j){
        while (i<j && a[j].value <= holeValue) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].value >= holeValue) i++;
        if (i<j){
            a[j]= a[i];
            j--;
        }
    }

    a[i].value = holeValue;
    a[i].id = holeId;
    angleQsort(l, i - 1, a);
    angleQsort(i + 1, r, a);
}

void pointQsort(int l, int r, POINT *a){
    if (l>=r)
        return;

    int i=l, j=r;
    POINT holeValue = a[i];
    while (i<j){
        while (i<j && a[j].x >= holeValue.x) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].x <= holeValue.x) i++;
        if (i<j){
            a[j] = a[i];
            j--;
        }
    }

    a[i] = holeValue;
    pointQsort(l  , i-1, a);
    pointQsort(i+1, r, a);
}

/*
 * rep0Cost1 = lambda0 * cost[0] + lambda1 * cost[1] + lambda2 * cost[2]
 * rep1Cost1 = lambda0 * cost[0] + lambda1 * cost[1] + lambda2 * cost[2]
 * ...
 * rep0Cost2 = lambda0' * cost[0] + lambda1' * cost[1] + lambda2' * cost[2]
 * rep1Cost2 = lambda0' * cost[0] + lambda1' * cost[1] + lambda2' * cost[2]
 * ...
 * rep0TotalCost  = rep0Cost1 + rep0Cost2 + ...
 */
void sortRepByLambda(myRep &rep){
    /*
     * apply for variables
     */
    int repSize = static_cast<int>(rep.size());
    EnergyValueStruct energyValue[nRep + 100];
    double **costInUnit;                //  costInUnit[1][2] : the third energy value of the second rep
    costInUnit = new double *[repSize];
    for (int i = 0; i < repSize; i++){
        costInUnit[i] =  new double [objectiveNumber];
    }

    getCostInUnit(rep, costInUnit);

    double *totalCost;
    totalCost = new double[repSize];        //  the total of new cost(add every energy value) for every rep
    memset(totalCost, 0, sizeof(totalCost));
    /*
     * 0 <= lambda <=1
     *   lamda[0][j] + lambda[1][j] + lamda[2][j] = 1
     */
    Matrix<double> lambda = randFixedSum(objectiveNumber, lambdaLoopTimes, 1.0, 0.0, 1.0);
    for (int j = 0; j < lambdaLoopTimes; j++){
        for (int i = 0; i < repSize; i++){
            double thisTimeValue = 0;
            for (int k = 0; k < objectiveNumber; k++){
                thisTimeValue += lambda[k][j] * costInUnit[i][k];
            }
            totalCost[i] += thisTimeValue;
        }
    }

    for (int i = 0; i < repSize; i++){
        cout << "i=" << i << " value=" << energyValue[i].value << endl;
        energyValue[i].value = totalCost[i] / lambdaLoopTimes;
        energyValue[i].id = i;
    }

    energyQsort(0, repSize - 1, energyValue);

    releaseSpace(costInUnit, repSize);
    releaseSpace(totalCost);


    cout << "sortAns :" << endl;
    for (int j = 0; j < repSize; j++){
        sortAns[j] = energyValue[j].id;
        cout << "j : " << j << " sortAns : " << sortAns[j] << endl;
    }
}

/*
 * from small to large
 * 1 2 3 4 5
 */
void energyQsort(int l, int r, EnergyValueStruct *a){
    if (l>=r)
        return;

    int i=l, j=r;
    double holeValue=a[i].value;
    int holeId = a[i].id;
    while (i<j){
        while (i<j && a[j].value >= holeValue) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].value <= holeValue) i++;
        if (i<j){
            a[j]= a[i];
            j--;
        }
    }

    a[i].value = holeValue;
    a[i].id = holeId;
    energyQsort(l, i - 1, a);
    energyQsort(i + 1, r, a);
}

void getCostInUnit(myRep &rep, double **costInUnit){
    int repSize = static_cast<int>(rep.size());
    /*
     * costInUnit = (x-min) / (max-min)
     */
    /*
     * get maxCost and minCost for 0~objectiveNumber
     */
    double maxCost[objectiveNumber], minCost[objectiveNumber];      // minCost[0] : the min value of the first energy function
    for (int i = 0; i < objectiveNumber; i++){
        maxCost[i] = -9999999999;
        minCost[i] = std::numeric_limits<double>::max();
    }
    myRep::iterator it = rep.begin();
    while (it != rep.end()) {
        for (int i = 0; i< objectiveNumber; i++){
            if (it->Cost[i] > maxCost[i])
                maxCost[i] = it->Cost[i];
            if (it->Cost[i] < minCost[i])
                minCost[i] = it->Cost[i];
        }
        it++;
    }
    /*
     * calculate newCost for every rep
     */
    it = rep.begin();
    int i = 0;
    while (it != rep.end()){
        for (int j = 0; j < objectiveNumber; j++)
            costInUnit[i][j] = (it->Cost[j] - minCost[j]) / (maxCost[j] - minCost[j]);
        i++;
        it++;
    }
}

template <typename _XD>
Matrix<double> randFixedSum(int n, int m, _XD S, _XD a, _XD b){
    //  n*a <= s <= n*b
    // Rescale to a unit cube: 0 <= x(i) <= 1
    S = S * 1.0;
    a = a * 1.0;
    b = b * 1.0;
    S = 1.0 * (S - n * a) / (b - a);
    /*
     * Construct the transition probability table, t.
     * t(i,j) will be utilized only in the region where j <= i + 1
     */
    int k = max( min(int(S), n-1), 0);     //   must have 0 <= k <= n-1
    S = max(min(S, 1.0*(k+1)), 1.0*k);               //   must have k <= s <= k+1

    Matrix<double> s1 = S - initArray(k, -1, k-n+1);
    Matrix<double> s2 = initArray(k+n, -1, k+1) - S;
    Matrix<double> w = Matrix<double>(n, n+1, 0);
    w[0][1] = realMax;                      // scale for full 'double' range
    Matrix<double> t = Matrix<double>(n-1, n, 0);

    for (int i = 1; i < n; i++){
        Matrix<double> tmp1 = mulCorNum(getArrayFromMatrix(w, i-1, 1, i+1+1),
                                        getArrayFromMatrix(s1, 0, 0, i+1)) / (i+1);
        Matrix<double> tmp2 = mulCorNum(getArrayFromMatrix(w, i-1, 0, i+1),
                                        getArrayFromMatrix(s2, 0, n-i-1, n-1+1)) / (i+1);
        w.matRowChange(i, 1, i+1+1, tmp1+tmp2);
        Matrix<double> tmp3 = getArrayFromMatrix(w, i, 1, i+1+1) + tiny;  // In case tmp1 & tmp2 are both 0
        Matrix<int> tmp4 = getArrayFromMatrix(s2, 0, n-i-1, n-1+1) >
                           getArrayFromMatrix(s1, 0, 0, i+1);                //  then t is 0 on left & 1 on right

        Matrix<double> tmpAns = mulCorNum(divCorNum(tmp2, tmp3), tmp4) +
                                mulCorNum((1 - divCorNum(tmp1, tmp3)), NOT(tmp4));
        t.matRowChange(i-1, 0, i+1, tmpAns);

    }
    Matrix<double> x(n, m, 0);
    Matrix<double> rt = randMatrix(n-1, m);     // For random selection of simplex type
    Matrix<double> rs = randMatrix(n-1, m);     // For random location withn simplex
    //rt = inputMatrix(n-1, m, 1.0);
    //rs = inputMatrix(n-1, m, 1.0);
    Matrix<double> s = repmat(S, 1, m);
    Matrix<int> j = repmat(k, 1, m);          //  For indexing in the t table
    Matrix<double> sm = Matrix<double>(1, m, 0.0);
    Matrix<double> pr = Matrix<double>(1, m, 1.0);          //  start with sum zero & product 1
    for (int i = n-2; i >= 0; i--){             // work backwards in the t table
        Matrix<int> e = getArrayFromMatrix(rt, n-i-2, 0, m) <= matRowBracket(t, i, j);  // Use rt to choose a transition
        Matrix<double> sx = matRoot(getArrayFromMatrix(rs, n-i-2, 0, m), 1.0/(i+1));        // Use rs to compute next simplex coord
        sm = sm + mulCorNum( mulCorNum(1-sx, pr), s) / (1.0*i+2);                       // Update sum
        pr = mulCorNum(sx, pr);                             // Update product
        x.matRowChange(n-i-2, 0, m, sm + mulCorNum(pr,e));  // Calculate x using simplex coords
        s = s - e;      j = j - e;                          // Transition adjustment
    }
    x.matRowChange(n-1, 0, m, sm + mulCorNum(pr, s));       // Compute the last x

    //  Randomly permute the order in the columns of x and rescale
    Matrix<double> rp = randMatrix(n, m);                   // Use rp to carry out a matrix 'randperm'
    //rp = inputMatrix(n, m, 1.0);
    Matrix<int> p = matrixColSort(rp);                      // The values placed in ig are ignored
    /*
     * x = (b-a)*x(p+repmat([0:n:n*(m-1)],n,1))+a   % Permute & rescale x
     * [0:n:n*(m-1)] , n = 3, m = 5  : [0 3 6 9 12]
     * repmat([0:n:n*(m-1)], n 1) : [0 3 6 9 12; 0 3 6 9 12; 0 3 6 9 12]  -- n*m
     */
    //cout << "Bracket:" << endl <<  p + repmat(initArray(0, n, n * (m-1)), n, 1) << endl;
    x = (b - a) * matBracket(x, p + repmat(initArray(0, n, n * (m-1)), n, 1)) + a;
    //cout << "x:" << endl << x << endl;
    return x;
}

/*
    * Sum of two Matrix
    */
template <typename _MD>
Matrix<_MD> operator+(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        throw std::invalid_argument("different matrix's size");
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<_MD> c(rowSize, colSize, 0);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++) {
            c[i][j] = a[i][j] + b[i][j];
        }
    return c;
}

template <typename _MD>
Matrix<_MD> operator+(const Matrix<_MD> &a, const _MD &b){              //  c[i][j] = a[i][j] + b;
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<_MD> c(rowSize, colSize, 0);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++) {
            c[i][j] = a[i][j] + b;
        }
    return c;
}

template <typename _MD>
Matrix<_MD> operator-(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        throw std::invalid_argument("different matrix's size");
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<_MD> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = a[i][j] - b[i][j];
    return c;
}

Matrix<double> operator-(const Matrix<double> &a, const Matrix<int> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        throw std::invalid_argument("different matrix's size");
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<double> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = a[i][j] - b[i][j];
    return c;
}

Matrix<double> operator-(const Matrix<int> &a, const Matrix<double> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        throw std::invalid_argument("different matrix's size");
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<double> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = a[i][j] - b[i][j];
    return c;
}

Matrix<double> operator-(const Matrix<int> &y, const double d){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<double> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = y[i][j] - d;
    return c;
}

Matrix<double> operator-(const Matrix<double> &y, const int d){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<double> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = y[i][j] - d;
    return c;
}

Matrix<double> operator-(const double d, const Matrix<int> &y){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<double> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = d - y[i][j];
    return c;
}

Matrix<double> operator-(const int d, const Matrix<double> &y){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<double> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = d - y[i][j];
    return c;
}

template <typename _MD>
Matrix<_MD> operator-(const Matrix<_MD> &y, const _MD d){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<_MD> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = y[i][j] - d;
    return c;
}

template <typename _MD>
Matrix<_MD> operator-(const _MD d, const Matrix<_MD> &y){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<_MD> c(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            c[i][j] = d - y[i][j];
    return c;
}

template <typename _MD>
bool operator==(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        return false;
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            if (a[i][j] != b[i][j])
                return false;
    return true;
}

template <typename _MD>
Matrix<_MD> operator-(const Matrix<_MD> &a){
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = -a[i][j];
    return result;
}

/*template <typename _MD>
Matrix<_MD> operator-(Matrix<_MD> &&a){
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            a[i][j] = -a[i][j];
    return a;
}*/
/*
 * multiplication of two matrix
 */
template <typename _MD>
Matrix<_MD> operator*(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.ColSize() != b.RowSize()){
        throw std::invalid_argument("can't multiply");
    }
    Matrix<_MD> result(a.RowSize(), b.ColSize(), 0);        // ***** forget 0
    size_t rowSize = a.RowSize(), colSize = b.ColSize(), pubSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            for (size_t k = 0; k < pubSize; k++)
                result[i][j] += a[i][k] * b[k][j];
    return result;
}
/*
 * matrix multiply a number
 */
template <typename _MD>
Matrix<_MD> operator*(const Matrix<_MD> &a, const _MD &b){
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = b * a[i][j];
    return result;
}

template <typename _MD>
Matrix<_MD> operator*(const _MD &b, const Matrix<_MD> &a){
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = b * a[i][j];
    return result;
}

template <typename _MD>
Matrix<_MD> operator/(const Matrix<_MD> &a, const _MD &b){
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] / b;
    return result;
}

Matrix<double> operator/(const Matrix<double> &a, const int &b){
    Matrix<double> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++) {
            result[i][j] = a[i][j] / b;
        }
    return result;
}

Matrix<double> operator/(const Matrix<int> &a, const double &b){
    Matrix<double> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] / b;
    return result;
}

/*
 * if a[i][j] > b[i][j], c[i][j] = 1, or c[i][j] = 0
 */
template <typename _MD>
Matrix<int> operator>(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        throw std::invalid_argument("different matrix's size");
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<int> c(rowSize, colSize, 0);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++) {
            c[i][j] = a[i][j] > b[i][j] ? 1 : 0;
        }
    return c;
}

/*
 * if a[i][j] <= b[i][j], c[i][j] = 1, or c[i][j] = 0
 */
template <typename _MD>
Matrix<int> operator<=(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()){
        throw std::invalid_argument("different matrix's size");
    }
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    Matrix<int> c(rowSize, colSize, 0);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++) {
            c[i][j] = a[i][j] <= b[i][j] ? 1 : 0;
        }
    return c;
}

/*   .*   规模相同的矩阵对应位置数相乘    mulCorNum(M,M)  //  multiply corresponding number
     ./  				             divCorNum(M,M)  //  two matrix with the same size
     .^				                 powCorNum(M,D)  */
template <typename _MD>
Matrix<_MD> mulCorNum(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.ColSize() != b.ColSize() || a.RowSize() != b.RowSize()){
        throw std::invalid_argument("can't multiply");
    }
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] * b[i][j];
    return result;
}

Matrix<double> mulCorNum(const Matrix<double> &a, const Matrix<int> &b){
    if (a.ColSize() != b.ColSize() || a.RowSize() != b.RowSize()){
        throw std::invalid_argument("can't multiply");
    }
    Matrix<double> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] * b[i][j];
    return result;
}

Matrix<double> mulCorNum(const Matrix<int> &a, const Matrix<double> &b){
    if (a.ColSize() != b.ColSize() || a.RowSize() != b.RowSize()){
        throw std::invalid_argument("can't multiply");
    }
    Matrix<double> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] * b[i][j];
    return result;
}

template <typename _MD>
Matrix<_MD> divCorNum(const Matrix<_MD> &a, const Matrix<_MD> &b){
    if (a.ColSize() != b.ColSize() || a.RowSize() != b.RowSize()){
        throw std::invalid_argument("can't divide");
    }
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] / b[i][j];
    return result;
}

Matrix<double> divCorNum(const Matrix<double> &a, const Matrix<int> &b){
    if (a.ColSize() != b.ColSize() || a.RowSize() != b.RowSize()){
        throw std::invalid_argument("can't divide");
    }
    Matrix<double> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] / b[i][j];
    return result;
}

Matrix<double> divCorNum(const Matrix<int> &a, const Matrix<double> &b){
    if (a.ColSize() != b.ColSize() || a.RowSize() != b.RowSize()){
        throw std::invalid_argument("can't divide");
    }
    Matrix<double> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = a[i][j] / b[i][j];
    return result;
}

template <typename _MD>
Matrix<_MD> powCorNum(const Matrix<_MD> &a, const int &p){
    Matrix<_MD> result(a.RowSize(), a.ColSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = matPow(a[i][j], p);
    return result;
}

/*
 * y = [1,2,3;4,5,6;7,8,9];
 * y(1,1:3) = [5,6];
 */
template <typename _MD>
Matrix<_MD> getArrayFromMatrix(const Matrix<_MD> &y, const int refRow, const size_t fromCol, const size_t endCol){
    size_t len = static_cast<size_t>(endCol - fromCol);
    Matrix<_MD> result(1, len);
    for (size_t j = fromCol; j < endCol; j++){
        result[0][j-fromCol] = y[refRow][j];
    }
    return result;
}

Matrix<int> initArray(const int &startNum, const int &delta, const int &endNum) {     //  realize [3:2：10] in matlab , need to be improved
    size_t len = 0;
    int k = startNum;
    while ((delta > 0 && k <= endNum) || (delta < 0 && k >= endNum)){
        len++;
        k += delta;
    }
    k = startNum;
    Matrix<int> result(1, len);
    for (size_t i = 0; i < len; i++){
        result[0][i] = k;
        k += delta;
    }
    return result;
}

template <typename _MD>         // repmat(MATRIX, n, m) /  repmat(double, n, m)
Matrix<_MD> repmat(const Matrix<_MD> &y, int n, int m){  //  copy Matrix for n times in row and m times in column
    size_t resRowSize = y.RowSize() * n, resColSize = y.ColSize() * m;
    size_t yRowSize = y.RowSize(), yColSize = y.ColSize();
    Matrix<_MD> result(resRowSize, resColSize);
    for (size_t yi = 0; yi < yRowSize; yi++){
        for (size_t yj = 0; yj < yColSize; yj++){
            for (size_t i = yi; i < resRowSize; i += yRowSize){
                for (size_t j = yj; j < resColSize; j += yColSize){
                    result[i][j] = y[yi][yj];
                }
            }
        }
    }
    return result;
}

Matrix<double> repmat(const double d, int n, int m){
    return Matrix<double>(n, m, d);
}

Matrix<int> repmat(const int d, int n, int m){
    return Matrix<int>(n, m, d);
}

template <typename _MD>
Matrix<_MD> transpose(const Matrix<_MD> &a){
    Matrix<_MD> result(a.ColSize(), a.RowSize());
    size_t rowSize = a.RowSize(), colSize = a.ColSize();
    for (size_t i = 0; i < colSize; i++)
        for (size_t j = 0; j < rowSize; j++)
            result[i][j] = a[j][i];
    return result;
}

//  unit MATRIX
template <typename _MD>
Matrix<_MD> I(const size_t &n){
    Matrix<_MD> result(n, n, 0);
    for (size_t i = 0; i < n; i++)
        result[i][i] = static_cast<_MD>(1);
    return result;
}

template <typename _MD>
Matrix<_MD> matPow(Matrix<_MD> mat, int b){         //  a^5
    if (mat.RowSize() != mat.ColSize()){
        throw std::invalid_argument("can't get power");
    }
    Matrix<_MD> result(I<_MD>(mat.RowSize()));
    while ( b > 0 ){
        if (b & 1)
            result = result*mat;
        mat = mat * mat;
        b >>= 1;
    }
    return result;
}

/*
 * Matrix mat;
 * mat^(1/2)
 */
template <typename _MD>
Matrix<_MD> matRoot(const Matrix<_MD> &mat, double p){
    size_t rowSize = mat.RowSize(), colSize = mat.ColSize();
    Matrix<_MD> result(rowSize, colSize);
    for (size_t i = 0; i < rowSize; i++)
        for (size_t j = 0; j < colSize; j++)
            result[i][j] = pow(mat[i][j], p);
    return result;
}


/*  matrixSort : sort by column
 *  [ig,p] = sort(rp);    p : the same size with rp, p[i][j] : p[i][j]是rp矩阵j列的第几个
 *  return p
 */
template <typename _MD>
Matrix<int> matrixColSort(const Matrix<_MD> y){
    Matrix<int> result(y.RowSize(), y.ColSize());
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    for (size_t j = 0; j < colSize; j++){
        commonStruct<_MD> a[rowSize];
        for (size_t i = 0; i < rowSize; i++){
            a[i].id = static_cast<int>(i);
            a[i].value = y[i][j];
        }
        colSort(0, int(rowSize-1), a);
        for (size_t i = 0; i < rowSize; i++){
            result[i][j] = a[i].id;
        }
    }
    return result;
}

template <typename _MD>
void colSort(int l, int r, commonStruct<_MD> *a){        // small first
    if (l>=r)
        return;

    int i=l, j=r;
    double holeValue=a[i].value;
    int holeId = a[i].id;
    while (i<j){
        while (i<j && a[j].value >= holeValue) j--;
        if (i<j){
            a[i] = a[j];
            i++;
        }
        while (i<j && a[i].value <= holeValue) i++;
        if (i<j){
            a[j]= a[i];
            j--;
        }
    }

    a[i].value = holeValue;
    a[i].id = holeId;
    colSort(l, i - 1, a);
    colSort(i + 1, r, a);
}

/*
 * MATRIX x, MATRIX y , z = x(y)
 * the same size with y, z[i][j] is the y[i][j](th) number in x Matrix(column first)
 */
template <typename _MD>
Matrix<_MD> matBracket(const Matrix<_MD> &x, const Matrix<int> &y){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<_MD> result(y.RowSize(), y.ColSize());
    for (size_t i = 0; i < rowSize; i++){
        for (size_t j = 0; j < colSize; j++){
            result[i][j] = getKthNumber(x, y[i][j]);
        }
    }
    return result;
}

/*
 * t = [1,2,3; 4,5,6];  j = [2,2,2];
 * t(0,j) = [3,3,3];
 * t(1,j) = [6,6,6]
 * this function is to obtain t(i,j); when j is a matrix
 */
template <typename _MD>
Matrix<_MD> matRowBracket(const Matrix<_MD> &x, const int &i, const Matrix<int> &y){
    size_t colSize = y.ColSize();
    Matrix<_MD> xRow = getArrayFromMatrix(x, i, 0, x.ColSize());
    Matrix<_MD> result(1,colSize, 0);
    for (size_t j = 0; j < colSize; j++){
        result[0][j] = getKthNumber(xRow, y[0][j]);
    }
    return result;
}

template <typename _MD>
_MD getKthNumber(const Matrix<_MD> &mat, int k){
    return mat[k - k/mat.RowSize() * mat.RowSize()][k/mat.RowSize()];
}

Matrix<double> randMatrix(const int n, const int m){
    Matrix<double> ans(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            ans[i][j] = rand() / double(RAND_MAX);
    return ans;
}

/*
 * change 1 to 0 / 0 to 1
 */
template <typename _MD>
Matrix<_MD> NOT(const Matrix<_MD> &y){
    size_t rowSize = y.RowSize(), colSize = y.ColSize();
    Matrix<_MD> result(y.RowSize(), y.ColSize());
    for (size_t i = 0; i < rowSize; i++){
        for (size_t j = 0; j < colSize; j++){
            result[i][j] = (y[i][j] + 1) % 2;
        }
    }
    return result;
}

/*
 * realize w(i,2:i+1) = tmp1 + tmp2 in matlab;
 * data[row][startCol:endCol-1] = x[0][~]
 * x must be an array Matrix
 */
template <typename _MD>
Matrix<_MD> Matrix<_MD>::matRowChange(const int &row, const int &startCol,
                                      const int &endCol, const Matrix<_MD> &x) {
    for (int j = startCol; j < endCol; j++){
        data[row][j] = x[0][j-startCol];
    }
    return *this;
}

template <typename _MD>
Matrix<_MD> inputMatrix(const size_t &n, const size_t &m, const _MD &type){
    Matrix<_MD> result(n, m);
    for (size_t i = 0 ; i < n; i++)
        for (size_t j = 0; j < m; j++)
            cin >> result[i][j];
    return result;
}

template <typename _MD>
std::ostream &operator<<(std::ostream &os, const Matrix<_MD> &mat){
    std::ostream::fmtflags oldFlags = os.flags();
    os.precision(8);
    os.setf(std::ios::fixed | std::ios::right);

    size_t rowSize = mat.RowSize(), colSize = mat.ColSize();
    os << std::endl;
    for (int i = 0; i < rowSize; i++){
        for (int j = 0; j < colSize; j++){
            os << std::setw(15) << mat[i][j] << " ";
        }
        os << std::endl;
    }

    os.flags(oldFlags);
    return os;
}