#include <iostream>
#include <vector>

#define ALMOST_ZERO 0.000005

using namespace std;

class Simplex {

public:
    int rows, cols;

    int colsOriginalA;
    int rowsOriginalA;

    std::vector<std::vector<double>> identityCertificate;
    std::vector<double> Certificate;

    std::vector<double> solution;

    double maximum = 0;

    int pivotColumToUseIlimited = -1;

    bool plIsIlimited;

    std::vector<std::vector<double>> A;
    std::vector<double> B;
    std::vector<double> C;

public:
    Simplex(std::vector<std::vector<double>> a, std::vector<double> b, std::vector<double> c, std::vector<std::vector<double>> certificate) {
        
        //inicializo valores
        maximum = 0;
        plIsIlimited = false;
        rows = a.size();
        cols = a[0].size();
        
        identityCertificate.resize(b.size(), vector<double>(b.size(), 0.0));
        Certificate.resize(b.size(), 0.0);

        B.resize(b.size());
        C.resize(c.size());
        A.resize(rows, vector<double>(cols, 0));

        for (int i = 0; i < rows; i++) { 
            for (int j = 0; j < rows; j++) {
                identityCertificate[i][j] = certificate[i][j];
            }
        }        

        for (int i = 0; i < rows; i++) { 
            for (int j = 0; j < cols; j++) {
                A[i][j] = a[i][j];
            }
        }

        for (int i = 0; i < b.size(); i++) { 
            B[i] = b[i];
        }

        for (int i = 0; i < c.size(); i++) { 
            C[i] = c[i];
        }

    }

    void prepareAux() { //preparo aux -> coloco na forma canonica

        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) { 
                C[j] -= A[i][j]; //zero variaveis com custo 1 representando identidade no vetor C

                if (j==0) {
                    Certificate[i] -= identityCertificate[i][i]; //atualizo certificado também
                    maximum -= B[i]; //atualizo otimo
                }
            }
        }        
        return;
    }
    
    bool simplexAlgorithmCalculataion() { //se for otimo -> paro. se não, faço pivoteamento caso não for ilimitada

        if (checkIfIsOptimal()) {
            return true;
        }

        vector<int> pivotLocation = findPivotLocation();

        if (plIsIlimited == true) {
            pivotColumToUseIlimited = pivotLocation[0];
            return true;
        }

        makeMatrixElimination(pivotLocation[1], pivotLocation[0]);

        return false;
    }

    bool checkIfIsOptimal() { //checo se vetor C é otimo
 
        for (int i = 0; i < C.size(); i++) {
            double value = C[i];

            if (value < 0) { //se tem algum valor < 0 então não é otimo
                return false;
            }
        }

        return true;

    }

    void makeMatrixElimination(int pivotRow, int pivotColumn) { //faço pivoteamento

        double pivotRowVals[cols]; //linha do pivo        
        double pivotRowValsCertificate[cols]; //coluna do pivo de certificado
        double pivotValue = A[pivotRow][pivotColumn]; 
        double rowNew[cols]; //linha do pivo atualizada
        double pivotColVals[rows]; //coluna do pivo
        double rowNewCertificate[rows]; // linha do pivo no certificado atualizada    

        //encontro novo maxio relativo ao pivo
        maximum = maximum - (C[pivotColumn] * (B[pivotRow] / pivotValue)); 

        for (int i = 0; i < cols; i++) { 
            pivotRowValsCertificate[i] = identityCertificate[pivotRow][i]; 
        }

        for (int m = 0; m < cols; m++) { 
            pivotRowVals[m] = A[pivotRow][m];
        }

        for (int i = 0; i < rows; i++) {
            pivotColVals[i] = A[i][pivotColumn];
        }

        //atualizo a linha do pivo
        for (int f = 0; f < cols; f++) {
            rowNew[f] = pivotRowVals[f] / pivotValue;

            if (abs(rowNew[f]) < ALMOST_ZERO ) rowNew[f] = 0; //arredondo valores muito pequenos para 0
        }

        //atualizo linha do pivo no certificado
        for (int k = 0; k < rows; k++) {
            rowNewCertificate[k] = pivotRowValsCertificate[k] / (double) pivotValue;

            if (abs(rowNewCertificate[k]) < ALMOST_ZERO) rowNewCertificate[k] = 0;  //arredondo valores muito pequenos para 0
        }

        B[pivotRow] = B[pivotRow] / pivotValue;

        if (abs(B[pivotRow]) < ALMOST_ZERO) B[pivotRow] = 0;  //arredondo valores muito pequenos para 0


        //calculo novos A, B e identidade certificado com base no pivo
        for (int i = 0; i < rows; i++) {

            if (pivotRow != i) {    //ignoro essa linha pois já calculei a linha do pivo
                    double valueToMultiplyRows = pivotColVals[i];

                    B[i] = B[i] - (valueToMultiplyRows * B[pivotRow]);
                    int invertSignal = 1;

                    if (abs(B[i]) < ALMOST_ZERO) B[i] = 0;  //arredondo valores muito pequenos para 0

                for (int j = 0; j < cols; j++) {
        
                    A[i][j] =  A[i][j] - (valueToMultiplyRows * rowNew[j]);

                    if (abs(A[i][j]) < ALMOST_ZERO) A[i][j] = 0;  //arredondo valores muito pequenos para 0

                    if (i < B.size() && j < B.size()) {
                        identityCertificate[i][j] = identityCertificate[i][j] - ((double) valueToMultiplyRows * rowNewCertificate[j]);

                        if (abs(identityCertificate[i][j]) < ALMOST_ZERO) identityCertificate[i][j] = 0;  //arredondo valores muito pequenos para 0
                    }
                }
            }
        }

        double valueToMultiplyRows = C[pivotColumn];
        for (int i = 0; i < C.size(); i++) { //atualizo custos
            C[i] = C[i] - (valueToMultiplyRows * rowNew[i]);

            if (abs(C[i]) < ALMOST_ZERO) C[i] = 0;
        }

        //atualizo certificado
        for (int i = 0; i < rows; i++) {
            Certificate[i] = Certificate[i] - ( valueToMultiplyRows * rowNewCertificate[i]);
            if (abs(Certificate[i]) < ALMOST_ZERO) Certificate[i] = 0;
        }

        //atualizo linha de A que era do pivo
        for (int i = 0; i < cols; i++) {
            A[pivotRow][i] = rowNew[i];
        }

        //atualizo linha do certificado que era do pivo
        for (int i = 0; i < rows; i++) {
            identityCertificate[pivotRow][i] = rowNewCertificate[i];
        }

        return;
    }

    vector<int> findPivotLocation () {
        int pivotCol = findColumnPivot();
        int pivotRow = findRowPivot(pivotCol);
        vector<int> pivotLocation;
        pivotLocation.push_back(pivotCol);
        pivotLocation.push_back(pivotRow);

        return pivotLocation;
    }

    int findColumnPivot() { //encontro primeiro valor negativo de C -> evita ciclagem

        for (int i = 0; i < C.size(); i++) {
            if (C[i] < 0.0) {
                return i;
            }
        }
        return -1;
    }

    int findRowPivot(int pivotColumn) {
        double positiveValues[rows];
        std::vector<double> result(rows, 0);
        int negativeValueCount = 0;

        for (int i = 0; i < rows; i++) {
            if (A[i][pivotColumn] > 0) { //salvo valores positivos na coluna do pivo
                positiveValues[i] = A[i][pivotColumn];
            }
            else {
                positiveValues[i] = 0;
                negativeValueCount += 1;
            }
        }
        // se toda linha de A no pivo é negativa -> ilimitada
        if (negativeValueCount == rows) {
            plIsIlimited = true;
            return 0;
        }
        else {
            for (int i = 0; i < rows; i++) {
                double value = positiveValues[i];
                if (value > 0 && B[i] >= 0) {
                    result[i] = B[i] / value; // vejo resultado da divisao dos valores de A por B para decidir de qual será linha será o pivo
                } else {
                    result[i] = -1;
                }
            }
        }

        double minimum = 99999999;
        int pivotRow = 0;

        for (int i = 0; i < rows; i++) {
            if (result[i] >= 0) {
                if (result[i] < minimum) { //acho qual será o pivo
                    minimum = result[i];

                    pivotRow = i;
                }
            }
        }

        return pivotRow;
    }

    int CalculatePLAux () {
        
        prepareAux();

        bool result = false;
        while (result == false) {

            result = simplexAlgorithmCalculataion();
        }

        //na PL auxiliar, só preciso do máximo para saber se é viável ou inviável
        return maximum;

    }

     std::vector<double> CalculateSimplex() {

        std::vector<double> CertificateIlimited(colsOriginalA, 0);
        
        bool result = false;
        while (result == false) {

            result = simplexAlgorithmCalculataion();
        }
        // após fim do simplex, encontro solução para variáveis e certificado (otimo ou ilimitada)

        for (int i = 0; i < colsOriginalA; i++) { 
            int index = 0;
            int hasBasicValueEquals1 = 0;
            int quantityValues0 = 0;
            for (int j = 0; j < rows; j++) {
                if (A[j][i] >= -0.001 && A[j][i] <= 0.001 ) { //conto quantidade de valores 0 para saber se coluna foi toda pivoteada
                    quantityValues0 += 1;
                }
                else if (A[j][i] >= 0.99 && A[j][i] <= 1.001) { // se elemento de A é 1 -> candidato para ser base
                    index = j;
                    hasBasicValueEquals1 = 1;
                }
            }

            //confirmo se a coluna atual é basica (tem unico valor 1 e o restante é 0)
            if ((rows - 1 == quantityValues0) && (1 == hasBasicValueEquals1)) {
                solution.push_back(B[index]); //salvo solucao -> referente ao valor de B associado a essa coluna basica

                if (plIsIlimited == true){ 
                    if (abs(A[index][i]) > ALMOST_ZERO) {
                        CertificateIlimited[i] = (-1) * A[index][pivotColumToUseIlimited];
                    } else {
                         CertificateIlimited[i] = 0.0;
                    }
                }
            }
            else { //se não for basica, solução é 0
                solution.push_back(0.0);
            }
        }

        //variavel que causou a PL ser ilimitada possui valor no certificado igual a 1
        if (plIsIlimited == true && pivotColumToUseIlimited < colsOriginalA) {
            CertificateIlimited[pivotColumToUseIlimited] = 1.0;
        }

        return CertificateIlimited;
    }

    std::vector<std::vector<double>> returnIdentityCertificateAux () {
        return identityCertificate;
    }

    std::vector<std::vector<double>> returnAAux () {
        return A;
    }

    std::vector<double> returnBAux () {
        return B;
    }


};

int main() {

    int colSizeA = -1; 
    int rowSizeA = -1;

    cin >> rowSizeA;
    cin >> colSizeA;

    double C[rowSizeA + colSizeA]; //custos reduzidos -> tem tamanho = quantidade linhas problema + var folga

    for (int i = 0; i < colSizeA; i++) {
        cin >>  C[i];
        C[i] = (-1.0) * C[i]; //Já deixo custos reduzidos negativos para o simplex
    }
    

    for (int i = colSizeA; i < rowSizeA + colSizeA; i++) {
        C[i] = 0.0; //custo reduzido das variaveis de folga
    }

    double B[rowSizeA]; 
    double a[rowSizeA][colSizeA+rowSizeA]; // crio A considerando variaveis de folga

    bool twoPhasesSimplex = false;

    std::vector<std::vector<double>> identityCertificate(rowSizeA, std::vector<double>(rowSizeA, 0));

    for (int i = 0; i < rowSizeA; i++) {
        identityCertificate[i][i] = 1.0; //inicializo identidade
    }

    //recebo valores da entrada
    for (int i = 0; i < rowSizeA; i++) {
        for (int j = 0; j < colSizeA; j++) {
            cin >> a[i][j];
        }
        cin >> B[i]; 

        for (int j = colSizeA; j < colSizeA+rowSizeA; j++) { // inicio variaveis de folga em A
            a[i][j] = 0.0;
            if (i == j - colSizeA) {
                a[i][j] = 1.0;
            }
        }

        if (B[i] < 0) { //trato caso em que B[i] é negativo
            twoPhasesSimplex = true; //se B tem valor negativo, então preciso fazer duas fases
            B[i] = B[i] * (-1.0);   //caso B negativo, multiplico certificado e linha de A por -1
            identityCertificate[i][i] = identityCertificate[i][i] * (-1.0);
            for (int j = 0; j < colSizeA+rowSizeA; j++) {
                a[i][j] = a[i][j] * (-1.0);
            }
        }        
    }
    
    //variaveis para fazer simplex duas fases -> pl auxiliar
    int colSizeAAux = colSizeA+rowSizeA; 
    std::vector<std::vector<double>> aAux(rowSizeA, std::vector<double>(colSizeAAux+rowSizeA, 0));
    std::vector<double> bAux(rowSizeA, 0);
    std::vector<double> cAux(colSizeAAux+rowSizeA, 0);

    //variaveis para o simplex na segunda fase
    std::vector<std::vector<double>> aSimplex(rowSizeA, std::vector<double>(colSizeA+rowSizeA, 0));
    std::vector<double> bSimplex(rowSizeA, 0);
    

    if (twoPhasesSimplex == true) {

        for (int i = 0; i < rowSizeA; i++) { // inicializo vetor A para pl auxiliar
            for (int j = 0; j < colSizeAAux; j++) {
                aAux[i][j] = a[i][j];
            }
            for (int j = colSizeAAux; j < colSizeAAux+rowSizeA; j++) { //criando identidade da auxiliar
                if (i == (j - colSizeAAux)) {
                    aAux[i][j] = 1.0;
                } else {
                    aAux[i][j] = 0.0;
                }
            }
        }
        for (int i = 0; i < rowSizeA; i++) { // Baux é copia do B
            bAux[i] = B[i];
        }

        for (int i = 0; i < colSizeAAux+rowSizeA; i++) { //Na auxiliar o C é todo zero com excecao das novas var da auxiliar
            cAux[i] = 0.0;
            if (i >= colSizeAAux) {
                cAux[i] = 1.0;
            }
        }        
        
        Simplex simplexAux(aAux, bAux, cAux, identityCertificate);
        simplexAux.colsOriginalA = colSizeA;
        simplexAux.rowsOriginalA = rowSizeA;

        int result = simplexAux.CalculatePLAux(); //faço pl auxiliar e retorno se é inviavel ou não
 
        if (result < 0) { //PL é inviável
            cout << "inviavel" << endl;
            for (auto i : simplexAux.Certificate) {
                printf("%.7f ", i);
            }
            cout << endl;
            return 0;
        
        }

        //Caso pl auxiliar não for inviavel -> Reaproveito identidade, B e A obtidos na pl auxiliar para fazer a segunda fase do simplex
        identityCertificate = simplexAux.returnIdentityCertificateAux();
        
        bSimplex = simplexAux.returnBAux();

        std::vector<std::vector<double>> vectorAux(rowSizeA, std::vector<double>(colSizeAAux+rowSizeA, 0));

        vectorAux = simplexAux.returnAAux();

        for (int i = 0; i < rowSizeA; i++) { 
            for (int j = 0; j < colSizeA+rowSizeA; j++) {
                aSimplex[i][j] = vectorAux[i][j];
            }
        }


    }

    std::vector<double> cSimplex(colSizeA+rowSizeA, 0);

    for (int i = 0; i < colSizeA+rowSizeA; i++) { //inicializo vetor C para fazer o simplex
            cSimplex[i] = C[i];
    }

    if (twoPhasesSimplex == false) { // caso não tenha feito duas fases, inicializo A e B para simplex
        for (int i = 0; i < rowSizeA; i++) { 
            for (int j = 0; j < colSizeA+rowSizeA; j++) {
                aSimplex[i][j] = a[i][j];
            }
        }

        for (int i = 0; i < rowSizeA; i++) {
            bSimplex[i] = B[i];
        }
        
    } 

    Simplex simplex(aSimplex, bSimplex, cSimplex, identityCertificate);
    simplex.colsOriginalA = colSizeA;
    simplex.rowsOriginalA = rowSizeA;

    std::vector<double> certificateIlimited(colSizeA, 0);   

    certificateIlimited = simplex.CalculateSimplex(); //faço simplex e retorno certificado de ilimitada (se houver)

    if (simplex.plIsIlimited == 0) { //PL tem valor otimo
        cout << "otima" << endl;
        printf("%.7f ", simplex.maximum);

        cout << endl;
        
        for (auto i : simplex.solution) {
            printf("%.7f ", i);
        }
        cout << endl;

        for (auto i: simplex.Certificate) {
            printf("%.7f ", i);
        }
        cout << endl;
    } else { //PL é ilimitada

        cout << "ilimitada" << endl;
        
        for (auto i : simplex.solution) {
            printf("%.7f ", i);
        }
        
        cout << endl;

        for (auto i: certificateIlimited) {
            printf("%.7f ", i);
        }
        cout << endl;
    }   

    return 0;
}