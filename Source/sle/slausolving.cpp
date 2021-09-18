#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"

#include "timeinterval.h"
#include "interpolation.h"

#include "slausolving.h"

namespace SlauSolving
{

const size_t maxBMPSize = 512;

// пустой предобусловливатель (отсутствие)
struct SSCMPreconditioner_none: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
};
// диагональный предобусловливатель (не полный)
struct SSCMPreconditioner_diag: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    // saveBMP ничего не делает
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    std::vector<double> d; // диагональные элементы
};
// LLT предобусловливатель, симметричный строчно-столбцовый формат, портрет не меняется (не полный)
struct SSCMPreconditioner_SSCM_LLT_0: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    SSCMElements m; // элементы матрицы L, партрет такой же как и у исходной
};
// LDLT предобусловливатель, симметричный строчно-столбцовый формат, портрет не меняется (не полный)
struct SSCMPreconditioner_SSCM_LDLT_0: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    SSCMElements m; // элементы матрицы L, партрет такой же как и у исходной
};
// LLT предобусловливатель, симметричный профильный формат, портрет расширяется (полный)
struct SSCMPreconditioner_Profile_LLT: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    ProfilePortrait p; // профильный портрет матрицы полного предобусловливателя
    SSCMElements e; // элементы матрицы L
};
// LDLT предобусловливатель, симметричный профильный формат, портрет расширяется (полный)
struct SSCMPreconditioner_Profile_LDLT: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    ProfilePortrait p; // профильный портрет матрицы полного предобусловливателя
    SSCMElements e; // элементы матрицы L и диагонали
};

// LLT предобусловливатель, симметричный строчно-столбцовый формат, портрет расширяется (полный)
struct SSCMPreconditioner_SSCM_LLT: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    SSCMPortrait p_SSCM; // строчно-столбцовый портрет матрицы полного предобусловливателя(может получиться расширенным)
    SSCMElements m_SSCM;
    //std::vector<double> s0; // строка - временный массив
    std::vector<double> s; // строка - временный массив
};
// LDLT предобусловливатель, симметричный строчно-столбцовый формат, портрет расширяется (полный)
struct SSCMPreconditioner_SSCM_LDLT: public SSCMPreconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector &x)const final;
    virtual void InversedQmulVector(Vector &x)const final;
private:
    SSCMPortrait p_SSCM; // строчно-столбцовый портрет матрицы полного предобусловливателя(может получиться расширенным)
    SSCMElements m_SSCM;
    std::vector<double> s; // строка - временный массив
};

// прямой решатель (предобусловие должно быть полное)
struct SSCMSolver_direct: public SSCMSolver_base
{
    virtual void init(const size_t matrixSize) final;
    virtual void release() final;
    virtual void solve(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, const Vector &x0, const SolverParameters &ssp,
               Vector &x, double &residual, double &relativeResidual, int &iterations, double &time) final;
private:
    Vector r;
    Vector r0;
};
// итерационный решатель (предобусловие может быть не полное)
struct SSCMSolver_iterative_base: public SSCMSolver_base
{
    virtual void solve(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, const Vector &x0, const SolverParameters &ssp,
                       Vector &x, double &residual, double &relativeResidual, int &iterations, double &time) final;
protected:
    virtual void start(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, Vector &x) = 0;
    virtual void iter(const SSCM &matrix, const Vector &, const SSCMPreconditioner_base *p, Vector &x) = 0;
    Vector r;
    Vector r0;
};
// итерационный решатель LOS
struct SSCMSolver_LOS: public SSCMSolver_iterative_base
{
    virtual void init(const size_t matrixSize) final;
    virtual void release() final;
private:
    virtual void start(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, Vector &x) final;
    virtual void iter(const SSCM &matrix, const Vector &, const SSCMPreconditioner_base *p, Vector &x) final;
    Vector ss;
    Vector z;
    Vector w;
    Vector aa;
    Vector pp;
};
// итерационный решатель CGM
struct SSCMSolver_CGM: public SSCMSolver_iterative_base
{
    virtual void init(const size_t matrixSize) final;
    virtual void release() final;
private:
    virtual void start(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, Vector &x) final;
    virtual void iter(const SSCM &matrix, const Vector &, const SSCMPreconditioner_base *p, Vector &x) final;
    Vector ss;
    Vector z;
    Vector w;
    Vector aa;
    Vector pp;
};

ProfilePortrait::ProfilePortrait(const size_t newMatrixSize)
{
    init(newMatrixSize);
}
void ProfilePortrait::init(const size_t newMatrixSize)
{
    matrixSize = newMatrixSize;
    ind.resize(matrixSize + 1);
}
void ProfilePortrait::init_LLT(const SSCMPortrait &p0)
{
    matrixSize = p0.getMatrixSize();
    // инициализация портрета
    size_t elementsNumber = 0; // общее количество элементов матрицы с плотным профилем
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        if(p0.ind[i + 1] - p0.ind[i] != 0)  // строка не пустая
        {
            size_t j0 = p0.ai[p0.ind[i]];   // j0 - номер столбца первого элемента строки i
            size_t stringSize = i - j0;     // в строке i плотный профиль занимает (i - j0) элементов
            elementsNumber += stringSize;
        }
    }
    ind.resize(matrixSize + 1);
    // заполнение профильного портрета
    ind[0] = 0;
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        if(p0.ind[i + 1] - p0.ind[i] != 0)  // строка не пустая
        {
            size_t j0 = p0.ai[p0.ind[i]];   // j0 - номер столбца первого элемента строки i
            size_t stringSize = i - j0;     // в строке i плотный профиль занимает (i - j0) элементов
            ind[i + 1] = ind[i] + stringSize;
        }
        else
            ind[i + 1] = ind[i];    // строка пустая
    }
}
void ProfilePortrait::release()
{
    matrixSize = 0;
    ind.clear();
}
void ProfilePortrait::saveBMP(const char *fileName, size_t size) const
{
    if(size > maxBMPSize)
        size = maxBMPSize;
    #define POINT(x_, y_) bmp.c[(bmp.x-(y_)-1)*bmp.x + (x_)]
    #define SCALED_POINT(x_, y_, xmax, ymax) POINT(((x_)*bmp.x/xmax), ((y_)*bmp.y/ymax))
    Interpolation::BMP_INF bmp;
    bmp.x = size;
    bmp.y = size;
    bmp.c = new unsigned int[bmp.x*bmp.y];
    for (int i = 0; i < bmp.x; i++)
        for (int j = 0; j < bmp.y; j++)
            POINT(j, i) = Interpolation::RGB32(0, 255, 255, 255);
    // вне диагонали
    for (size_t i = 0; i < matrixSize; i++)           // i - номер строки
    {
        size_t stringSize = ind[i + 1] - ind[i];
        for (size_t j = i - stringSize; j <= i; j++)           // j - номер столбца
        {
            if(SCALED_POINT(j, i, matrixSize, matrixSize) != 0)
            {
                SCALED_POINT(j, i, matrixSize, matrixSize) = 0;
                //SCALED_POINT(i, j, matrixSize, matrixSize) = 0;
            }
        }
    }
    // сохранение картинки
    Interpolation::bmp_save(fileName, bmp);
    delete[] bmp.c;
}
void ProfilePortrait::saveProperties(const char *fileName)const
{
    FILE *f = fopen(fileName, "w");
    fprintf(f, "Size = %dx%d = %.0lf\n", (int)getMatrixSize(), (int)getMatrixSize(), SQR((double)getMatrixSize()));
    fprintf(f, "Filled = %d, %.4lf percents\n", (int)getElementsNumber(), (double)getElementsNumber()/SQR((double)getMatrixSize())*100);
    fclose(f);
}
size_t ProfilePortrait::getMatrixSize()const
{
    return matrixSize;
}
size_t ProfilePortrait::getElementsNumber()const
{
    return ind[matrixSize];
}

SSCMPortrait::SSCMPortrait(const size_t newElementsNumber, const size_t newMatrixSize)
{
    init(newElementsNumber, newMatrixSize);
}
void SSCMPortrait::init(const size_t newElementsNumber, const size_t newMatrixSize)
{
    matrixSize = newMatrixSize;
    ind.resize(matrixSize + 1);
    elementsNumber = newElementsNumber;
    ai.resize(elementsNumber);
}

// #можно добавить список очевидно заполенных ячеек(номеров столбцов, упорядоченных)
// и параллельно идти по этому списку (а не по портрету исходной матрицы)
// но наверно будет грузиться кеш
void SSCMPortrait::init_LLT(const SSCMPortrait &p0)
{
TimeIntervals::timeInterval debug_t;
debug_t.begin();
    matrixSize = p0.matrixSize;
    // 1) построение портрета разложения p (символическое разложение)
    ind.resize(matrixSize + 1);
    ind[0] = 0;
    ai.clear();
    ai.reserve(p0.getElementsNumber());
int lastPerc = 0;
double lastTime = 0;
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        // нужно построить портрет i-й строки, то есть ai[k] при ind[i]<=k<ind[i+1] - должны стать равны номерам столбцов элементов
        // обнуление строки
        ind[i + 1] = ind[i];
        if(p0.ind[i + 1] - p0.ind[i] != 0)
        {
            size_t k0 = p0.ind[i];
            size_t j0 = p0.ai[k0];  // номер столбца ненулевого элемента в портрете исходной матрицы
            size_t j_start = j0;  // номер столбца первого ненулевого элемента текущей строки в портрете исходной матрицы
            // проход по плотному портрету текущей строки, проверка пустых ячеек
            for (size_t j = j_start; j < i; j++) // j - номер столбца
            {
                if(j0 == j || SSCMintersection(*this, i, j, ind[i + 1]))
                {
                    ai.push_back(j);
                    //ai[ind[i + 1]] = j;
                    ind[i + 1]++;
                }
                // параллельно перемещаемся по исходному портрету
                if(j >= j0 && k0 + 1 < p0.ind[i + 1])
                {
                    k0++;
                    j0 = p0.ai[k0];
                }
            }
        }
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
  //fprintf(stderr, "%d\n", (int)ai.max_size());
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)p0.ind[i + 1]/p0.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    elementsNumber = ind[matrixSize];
}

/*void SSCMPortrait::init_LLT(const SSCMPortrait &p0)
{
TimeIntervals::timeInterval debug_t;
debug_t.begin();
    matrixSize = p0.matrixSize;
    // 1) построение портрета разложения p (символическое разложение)
    ind.resize(matrixSize + 1);
    ind[0] = 0;
    ai.clear();
    ai.reserve(p0.getElementsNumber());
    // вспомогательные массивы
    std::vector<size_t> firstNotNull(matrixSize, 0);
    // firstNotNull[j] - номер первой строки (начиная с диагонали), в которой ячейка столбца j ненулевая
    // Если ячейка (i, j) ненулевая и firstNotNull[j] != 0, то ячейка (i, firstNotNull[j]) - заполненная
    // инициализация 0 - во всех столбцах ненулевые элементы под диагональю ещё не встретились
    std::vector<bool> stringPortrait(matrixSize, false); // портрет строки
    // stringPortrait[j - j_start] = true если ячейка (i, j) ненулевая
    // очистка плотного портрета текущей строки
int lastPerc = 0;
double lastTime = 0;
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        // нужно построить портрет i-й строки, то есть ai[k] при ind[i]<=k<ind[i+1] - должны стать равны номерам столбцов элементов
        // обнуление строки
        ind[i + 1] = ind[i];
        if(p0.ind[i + 1] - p0.ind[i] != 0)
        {
            size_t j_min = p0.ai[p0.ind[i]];  // номер столбца первого ненулевого элемента в портрете текущей строки исходной матрицы
            // проход по портрету строки исходной матрицы - заполнение очевидно ненулевых ячеек
            for(size_t k0 = p0.ind[i]; k0 < p0.ind[i + 1]; k0++)
            {
                size_t j0 = p0.ai[k0];  // номер столбца ненулевого элемента в портрете исходной матрицы
                // добавление в портрет ячейки (i, j0)
                stringPortrait[j0] = true;
                if(firstNotNull[j0] != 0)
                {
                    // добавление в портрет ячейки (i, firstNotNull[j0])
                    stringPortrait[firstNotNull[j0]] = true;
                }
                else
                {
                    // i - номер первой строки, в которой ячейка столбца j0 ненулевая
                    firstNotNull[j0] = i;
                }
            }
            size_t j_max = j_min; // последний ненулевой элемент текущей строки (на данный момент)
            // проход по плотному портрету текущей строки, проверка пустых ячеек
            for (size_t j = j_min; j < i; j++) // j - номер столбца
            {
                // interception(p, i, j) - проверка наличия пересечения портретов строк i и j
                //if(stringPortrait[j] == true || SSCMintersection(p, i, j))
                if(stringPortrait[j] == true || SSCMintersection_fast(j, *this, stringPortrait, j_min, j_max))
                {
                    ai.push_back(j);
                    j_max = j;
                    //ai[ind[i + 1]] = j;
                    ind[i + 1]++;
                    if(firstNotNull[j] == 0)
                    {
                        // i - номер первой строки, в которой ячейка столбца j ненулевая
                        firstNotNull[j] = i;
                    }
                    stringPortrait[j] = true;
                }
            }
            // очистка плотного портрета текущей строки
            for(size_t k = ind[i]; k < ind[i + 1]; k++)
            {
                // ai[k] - номер столбца ненулевого элемента в портрете текущей строки матрицы разложения
                stringPortrait[ai[k]] = false; // обнуление
            }
        }
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
  //fprintf(stderr, "%d\n", (int)ai.max_size());
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)p0.ind[i + 1]/p0.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    elementsNumber = ind[matrixSize];
}
*/
/*
size_t j_start = p0.ai[p0.ind[i]];  // номер столбца первого ненулевого элемента в портрете текущей строки исходной матрицы
// очистка плотного портрета текущей строки
for (size_t j = j_start; j < i; j++) // j - номер столбца
{
    stringPortrait[j - j_start] = false;
}
// проход по портрету строки исходной матрицы - заполнение очевидно ненулевых ячеек
for(size_t k0 = p0.ind[i]; k0 < p0.ind[i + 1]; k0++)
{
    size_t j0 = p0.ai[k0];  // номер столбца ненулевого элемента в портрете исходной матрицы
    // добавление в портрет ячейки (i, j0)
    stringPortrait[j0 - j_start] = true;
    if(firstNotNull[j0] != 0)
    {
        // добавление в портрет ячейки (i, firstNotNull[j0])
        stringPortrait[firstNotNull[j0] - j_start] = true;
    }
    else
    {
        // i - номер первой строки, в которой ячейка столбца j0 ненулевая
        firstNotNull[j0] = i;
    }
}
// проход по плотному портрету текущей строки, проверка пустых ячеек
for (size_t j = j_start; j < i; j++) // j - номер столбца
{
    // interception(p, i, j) - проверка наличия пересечения портретов строк i и j
    if(stringPortrait[j - j_start] == true || SSCMintersection(p, i, j))
    {
        ai.push_back(j);
        //ai[ind[i + 1]] = j;
        ind[i + 1]++;
        if(firstNotNull[j] == 0)
        {
            // i - номер первой строки, в которой ячейка столбца j ненулевая
            firstNotNull[j] = i;
        }
    }
}
*/
void SSCMPortrait::release()
{
    matrixSize = 0;
    elementsNumber = 0;
    ind.clear();
    ai.clear();
}
void SSCMPortrait::saveBMP(const char *fileName, size_t size) const
{
    if(size > maxBMPSize)
        size = maxBMPSize;
    #define POINT(x_, y_) bmp.c[(bmp.x-(y_)-1)*bmp.x + (x_)]
    #define SCALED_POINT(x_, y_, xmax, ymax) POINT(((x_)*bmp.x/xmax), ((y_)*bmp.y/ymax))
    Interpolation::BMP_INF bmp;
    bmp.x = size;
    bmp.y = size;
    bmp.c = new unsigned int[bmp.x*bmp.y];
    for (int i = 0; i < bmp.x; i++)
        for (int j = 0; j < bmp.y; j++)
            POINT(j, i) = Interpolation::RGB32(0, 255, 255, 255);
    // вне диагонали
    for (size_t i = 0; i < matrixSize; i++)           // i - номер строки
        for (size_t t = ind[i]; t < ind[i + 1]; t++)
        {                               // t - индекс элемента i-й строки
            size_t j = ai[t];              // j - номер столбца
                                        // A(i, j) = e.a[t]
            if(SCALED_POINT(j, i, matrixSize, matrixSize) != 0)
            {
                SCALED_POINT(j, i, matrixSize, matrixSize) = 0;
                //SCALED_POINT(i, j, matrixSize, matrixSize) = 0;
            }
        }
    // диагональ
    for (size_t i = 0; i < matrixSize; i++)
        SCALED_POINT(i, i, matrixSize, matrixSize) = 0;
    // сохранение картинки
    Interpolation::bmp_save(fileName, bmp);
    delete[] bmp.c;
}
void SSCMPortrait::saveProperties(const char *fileName)const
{
    FILE *f = fopen(fileName, "w");
    fprintf(f, "Size = %dx%d = %.0lf\n", (int)getMatrixSize(), (int)getMatrixSize(), SQR((double)getMatrixSize()));
    fprintf(f, "Filled = %d, %.4lf percents\n", (int)getElementsNumber(), (double)getElementsNumber()/SQR((double)getMatrixSize())*100);
    fclose(f);
}
size_t SSCMPortrait::getMatrixSize()const
{
    return matrixSize;
}
size_t SSCMPortrait::getElementsNumber()const
{
    return ind[matrixSize];
}
size_t SSCMPortrait::findUnsorted(const size_t i, const size_t j)const
{
    for(size_t count_a = ind[i]; count_a < ind[i + 1]; count_a++)
    {
        if(ai[count_a] == j)
            return count_a;
    }
    return -1;
}
size_t SSCMPortrait::findSorted(const size_t i, const size_t j)const
{
    size_t min = ind[i];
    size_t max = ind[i + 1] - 1;
    for(;;)
    {
        size_t count_a = (min + max) / 2;
        if(ai[count_a] == j)
            return count_a;
        if(ai[count_a] < j)
            min = count_a + 1;
        else
            max = count_a - 1;
        if(min > max)
            break;
    }
    return -1;
}

SSCMElements::SSCMElements(const size_t elementsNumber, const size_t matrixSize)
{
    init(elementsNumber, matrixSize);
}
void SSCMElements::init(const size_t elementsNumber, const size_t matrixSize)
{
    a.resize(elementsNumber);
    d.resize(matrixSize);
}
void SSCMElements::fill(const double value)
{
    for(size_t i = 0; i < a.size(); i++)
        a[i] = value;
    for(size_t i = 0; i < d.size(); i++)
        d[i] = value;
}
void SSCMElements::release()
{
    a.clear();
    d.clear();
}

void SSCM::init()
{
    p = new SSCMPortrait;
    e = new SSCMElements;
}
void SSCM::release()
{
    p->release();
    delete p;
    e->release();
    delete e;
}
size_t SSCM::getMatrixSize()const
{
    return p->matrixSize;
}
size_t SSCM::getElementsNumber()const
{
    return p->elementsNumber;
}

void SSCMcopy(const SSCMElements &e1, SSCMElements &e2)
{
    if(e2.a.size() < e1.a.size())
        e2.a.resize(e1.a.size());
    if(e2.d.size() < e1.d.size())
        e2.d.resize(e1.d.size());
    for (size_t i = 0; i < e1.a.size(); i++)
        e2.a[i] = e1.a[i];
    for (size_t i = 0; i < e1.d.size(); i++)
        e2.d[i] = e1.d[i];
}
void SSCMmulScalar(SSCM &matrix, const double c)
{
    SSCMElements &e = *matrix.e;
    for (size_t i = 0; i < e.a.size(); i++)
        e.a[i] *= c;
    for (size_t i = 0; i < e.d.size(); i++)
        e.d[i] *= c;
}
void SSCMmulVector(const SSCM &matrix, const Vector &x, Vector &y)
{
    SSCMPortrait &p = *matrix.p;
    SSCMElements &e = *matrix.e;
    for (size_t i = 0; i < p.matrixSize; i++)
        y[i] = 0;
    // вне диагонали
    for (size_t i = 0; i < p.matrixSize; i++)           // i - номер строки
        for (size_t k = p.ind[i]; k < p.ind[i + 1]; k++)
        {                             // k - индекс элемента i-й строки
            size_t j = p.ai[k];          // j - номер столбца
                                      // A(i, j) = e.a[k]
            y[i] += e.a[k] * x[j];    // нижний треугольник
            y[j] += e.a[k] * x[i];    // верхний треугольник
        }
    // диагональ
    for (size_t i = 0; i < p.matrixSize; i++)
        y[i] += e.d[i] * x[i];
}
void SSCMsolveResidual(const SSCM &matrix, const Vector &b, const Vector &x, Vector &r)
{
    SSCMmulVector(matrix, x, r);            // A*x -> r
    Vector1PlusCmulVector2(b, -1, r, r);    // b - r -> r
}
void SSCMaddBoundaryCondition1(SSCM &matrix, Vector &b, const std::vector<double>&u0, const std::vector<bool> &state)
{
    SSCMPortrait &p = *matrix.p;
    SSCMElements &e = *matrix.e;
    // вне диагонали
    for (size_t i = 0; i < p.matrixSize; i++)           // i - номер строки
        for (size_t k = p.ind[i]; k < p.ind[i + 1]; k++)
        {                             // k - индекс элемента i-й строки
            size_t j = p.ai[k];          // j - номер столбца
                                      // A(i, j) = e.a[k]
            if(state[i])
            {
                b[j] -= u0[i]*e.a[k];   //j - номер строки зеркального элемента к e.a[k]
                e.a[k] = 0;
            }
            if(state[j])
            {
                b[i] -= u0[j]*e.a[k];
                e.a[k] = 0;
            }
        }
    // диагональ
    for (size_t i = 0; i < p.matrixSize; i++)
    {
        if(state[i])
        {
            e.d[i] = 1;
            b[i] = u0[i];
        }
    }
}
void SSCMaddBoundaryCondition1(SSCM &matrix, Vector &b, const size_t i0, const double t)
{
    SSCMPortrait &p = *matrix.p;
    SSCMElements &e = *matrix.e;
    size_t matrixSize = p.matrixSize;
    e.d[i0] = 1;
    b[i0] = t;
    for (size_t j = p.ind[i0]; j < p.ind[i0 + 1]; j++)
    {
        b[p.ai[j]] -= t*e.a[j]; //ai[j] - номер строки зеркального элемента к a[j]
        e.a[j] = 0;
    }
    // строка i0 учла первое краевое условие
    // осталось вычесть ее из остальных строк чтобы сохранить симметричность матрицы
    for (size_t i = i0 + 1; i < matrixSize; i++)       // идем сверху вниз по строкам i
    {
        int count1 = p.ind[i];
        int count2 = p.ind[i + 1] - 1;
        // поиск p.ai[count] = i0 для count1 <= count <= count2
        for (;;)
        {
            int count = (count1 + count2) / 2;
            if (p.ai[count] < i0)
                count1 = count + 1;
            else
                if (p.ai[count] > i0)
                    count2 = count - 1;
                else
                {// найден
                    b[i] -= t*e.a[count];
                    e.a[count] = 0;
                    break;
                }
            if (count1 > count2) break;
        }
    }
}
void SSCM1addSimilarM2(SSCM &matrix1, SSCM &matrix2)
{
    SSCMElements &e1 = *matrix1.e;
    SSCMElements &e2 = *matrix2.e;
    for (size_t i = 0; i < e1.a.size(); i++)
        e1.a[i] += e2.a[i];
    for (size_t i = 0; i < e1.d.size(); i++)
        e1.d[i] += e2.d[i];
}
void SSCM1addEnclosedM2mulScalar(SSCM &matrix1, SSCM &matrix2, const double c)
{
    SSCMPortrait &p1 = *matrix1.p;
    SSCMElements &e1 = *matrix1.e;
    SSCMPortrait &p2 = *matrix2.p;
    SSCMElements &e2 = *matrix2.e;
    //FILE *out = fopen("_", "w");
    size_t matrixSize = p1.matrixSize;         // размеры матриц должны быть одинаковыми
    // вне диагонали
    for (size_t i = 0; i < matrixSize; i++)    // i - номер строки обоих матриц
    {
        size_t k1 = p1.ind[i];                 // k1 - индекс элемента i-й строки матрицы 1
        size_t k2 = p2.ind[i];                 // k2 - индекс элемента i-й строки матрицы 2
        for (; k2 < p2.ind[i + 1]; k2++)    // обходим элементы матрицы 2
        {
            size_t j2 = p2.ai[k2];             // j2 - номер столбца элемента e2.a[k2]
                                            // M2(i, j2) = e2.a[k2]
            // поиск элемента (i, j2) в матрице M1
            // он должен существовать ввиду вложенности портретов поэтому проверка k1 < p1.ind[i + 1] не нужна
            for(;;)
            {
                size_t j1 = p1.ai[k1];   // j1 - номер столбца элемента e1.a[k1]
                                      // M1(i, j1) = e1.a[k1]
                if(j1 == j2)
                {
                    // M1(i, j1) += M2(i, j2)
                    e1.a[k1] += e2.a[k2] * c;
                    k1++;
                    break;
                }
                k1++;
            }
        }
    }
    //fclose(out);
    // диагональ
    for (size_t i = 0; i < matrixSize; i++)
        e1.d[i] += e2.d[i] * c;
}

// копирование строки i из матрицы p2, e2 в матрицу p1, e1
// матрицы должны быть одного размера, p2 вложен в p1
// элементы которых нет в портрете p2 задаются нулями в матрице p1
void SSCMStrCopyToSSCMStr(const size_t i, const SSCMPortrait &p1, SSCMElements &e1, const SSCMPortrait &p2, const SSCMElements &e2)
{
    // вне диагонали
    size_t k1 = p1.ind[i];          // k1 - индекс элемента i-й строки матрицы 1
    size_t k1_max = p1.ind[i + 1];
    size_t k2 = p2.ind[i];          // k2 - индекс элемента i-й строки матрицы 2
    size_t k2_max = p2.ind[i + 1];
    while(k2 < k2_max) // обход элементов строки матрицы 2
    {
        size_t j2 = p2.ai[k2];      // j2 - номер столбца элемента e2.a[k2]
                                    // M2(i, j2) = e2.a[k2]
        // поиск элемента (i, j2) в матрице M1
        // он должен существовать ввиду вложенности портретов поэтому проверка k1 < p1.ind[i + 1] не нужна
        for(;;) // обход элементов строки матрицы 1
        {
            size_t j1 = p1.ai[k1];  // j1 - номер столбца элемента e1.a[k1]
                                    // M1(i, j1) = e1.a[k1]
            if(j1 == j2)
            {
                // M1(i, j1) = M2(i, j2)
                e1.a[k1] = e2.a[k2];
                k1++;
                k2++;
                break;
            }
            else
            {
                e1.a[k1] = 0;
                k1++;
            }
        }
    }
    // заполнение нулями остатка строки матрицы M1
    while(k1 < k1_max)  // обход элементов матрицы 1
    {
        e1.a[k1] = 0;
        k1++;
    }
    // диагональ
    e1.d[i] = e2.d[i];
}

// копирование строки i из SSCM матрицы p2, e2 в профильную матрицу p1, e1
// матрицы должны быть одного размера, p2 вложен в p1
// элементы которых нет в портрете p2 задаются нулями в матрице p1
void SSCM3x3StrCopyToProfileStr(const size_t i, const ProfilePortrait &p1, SSCMElements &e1, const SSCMPortrait &p2, const SSCMElements &e2)
{
    // копия матрицы
    if(p2.ind[i + 1] - p2.ind[i] != 0)  // строка не пустая
    {
        size_t k1_min = p1.ind[i];
        size_t k1_max = p1.ind[i + 1];
        size_t k2_min = p2.ind[i];
        size_t k2_max = p2.ind[i + 1];
        size_t j_min = p1.getFirstCol(i);// j_min - номер столбца первого ненулевого элемента строки i
        //size_t stringSize = i - j_min;     // в строке i плотный профиль занимает stringSize элементов
        for(size_t k1 = k1_min; k1 < k1_max; k1++)
            e1.a[k1] = 0;     // e1.a[k] = L(i, i-stringSize + k-p1.ind[i])
        for(size_t k2 = k2_min; k2 < k2_max; k2++)
        {
            size_t j = p2.ai[k2];           // e2.a[k2] = A(i, j)
            size_t k1 = k1_min + (j - j_min); // e1.a[k1] = L(i, j)
            e1.a[k1] = e2.a[k2];
        }
    }
    e1.d[i] = e2.d[i];
}

// возвращает true если строки i1 и i2 пересекаются (i1 ниже чем i2)
// странно, но работает чуть быстрее чем SSCMintersection_fast (т.к. там куча массивов грузит кэш)
bool SSCMintersection(const SSCMPortrait &p, const size_t i1, const size_t i2, const size_t k1_max)
{
    size_t k1 = k1_max;//p.ind[i1 + 1];
    size_t k1_min = p.ind[i1];
    size_t k2 = p.ind[i2 + 1];
    size_t k2_min = p.ind[i2];
    if(k1 == k1_min ||
       k2 == k2_min)
        return false;
    k1--;
    k2--;
    //size_t j1_min = p.ai[k1_min];
    size_t j1_max = p.ai[k1];
    size_t j2_min = p.ai[k2_min];
    size_t j2_max = p.ai[k2];
    if(j1_max < j2_min /*|| j2_max < j1_min*/)
        return false;
    size_t j1 = j1_max;  // номер столбца в строке i1
    size_t j2 = j2_max;  // номер столбца в строке i2
    for(;;)
    {
        if(j1 == j2)
        {
            return true;
        }
        else
        {
            if(j1 > j2)
            {
                for(;;)
                {
                    if(k1 == k1_min)
                        return false;
                    k1--;
                    j1 = p.ai[k1];
                    if(j1 <= j2)
                        break;
                }
            }
            else
            {
                for(;;)
                {
                    if(k2 == k2_min)
                        return false;
                    k2--;
                    j2 = p.ai[k2];
                    if(j1 >= j2)
                        break;
                }
            }
        }
    }
    return false;
    /*
    size_t k1 = p.ind[i1];
    size_t k1_max = p.ind[i1 + 1];
    size_t k2 = p.ind[i2];
    size_t k2_max = p.ind[i2 + 1];
    if(k1_max - k1 == 0 ||
       k2_max - k2 == 0)
        return false;
    for(;;)
    {
        size_t j1 = p.ai[k1];  // номер столбца в строке i1
        size_t j2 = p.ai[k2];  // номер столбца в строке i2
        if(j1 == j2)
        {
            return true;
        }
        else
        {
            if(j1 < j2)
            {
                k1++;
                if(k1 == k2_max)
                    break;
            }
            else
            {
                k2++;
                if(k2 == k2_max)
                    break;
            }
        }
    }
    return false;*/
}

// возвращает true если портрет строки i матрицы (с портретом p) пересекает строку с портретом s
/*bool SSCMintersection_fast(const size_t i, const SSCMPortrait &p, const std::vector<bool> &s, const size_t s_j_min, const size_t s_j_max)
{
    size_t k_min = p.ind[i];
    size_t k_max = p.ind[i + 1];
    if(k_min == k_max)
        return 0;
    size_t m_j_min = p.ai[k_min];
    size_t m_j_max = p.ai[k_max - 1];
    if(s_j_max < m_j_min || m_j_max < s_j_min)
        return false;
    size_t j_min = MAX(s_j_min, m_j_min);   // столбцы j < j_min нулевые для строки матрицы и строки вектора s
    size_t k = k_max - 1;
    for(;;) // обход элементов строки
    {
        size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
        if(s[j])
            return true;
        if(j <= j_min)
            return false;
        k--;
    }
}
*/

// вычисление скалярного произведения строки i матрицы (p, e) на плотный вектор s
// s_j_min - первый ненулевой столбец строки, соответствующей вектору s
double SSCMStrScalarMulVector(const size_t i, const SSCMPortrait &p, const SSCMElements &e, const std::vector<double> &s, const size_t s_j_min)
{
    // слева направо
    /*
    double E = 0;
    size_t k = p.ind[i]; // k - индекс элемента i-й строки
    size_t k_max = p.ind[i + 1];
    while(k < k_max) // обход элементов строки
    {
        if(s[p.ai[k]] != 0)
        E += e.a[k] * s[p.ai[k]];// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
        k++;
    }
    return E;
    */
    // справа налево
    double E = 0;
    size_t k_min = p.ind[i];
    size_t k_max = p.ind[i + 1];
    if(k_min == k_max)
        return 0;
    size_t m_j_min = p.ai[k_min];
    size_t j_min = MAX(s_j_min, m_j_min);   // столбцы j < j_min нулевые для строки матрицы и строки вектора s
    size_t k = k_max - 1;

    for(;;) // обход элементов строки
    {
        size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
        //if(s[j] != 0)
            E += e.a[k] * s[j];// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
        if(j <= j_min)
            break;
        k--;
    }

    // умножаем пачками по size
    /*
#define NNNN 32
    //const size_t size = 32;
    for(;;) // обход элементов строки
    {
        bool defaultCycle = false;
        if(k >= k_min + NNNN)
        {
            if(p.ai[k - NNNN] > j_min)
            {
                for(int i = 0; i < NNNN; i++)
                {
                    size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
                    E += e.a[k] * s[j];// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
                    k--;
                }
            }
            else
            {
                defaultCycle = true;
            }
        }
        else
        {
            defaultCycle = true;
        }

        defaultCycle = true;

        if(defaultCycle)
        {
            for(;;)
            {
                size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
                E += e.a[k] * s[j];// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
                if(j <= j_min)
                    goto end_SSCMStrScalarMulVector;
                k--;
            }
        }
    }
    end_SSCMStrScalarMulVector:;
    */

    return E;
}

// вычисление скалярного произведения строки i матрицы (p, e) на плотный вектор s, домножается дополнительно на e.d[]
// s_j_min - первый ненулевой столбец строки, соответствующей вектору s
double SSCMStrScalarMulVectoMulD(const size_t i, const SSCMPortrait &p, const SSCMElements &e, const std::vector<double> &s, const size_t s_j_min)
{
    // справа налево
    double E = 0;
    size_t k_min = p.ind[i];
    size_t k_max = p.ind[i + 1];
    if(k_min == k_max)
        return 0;
    size_t m_j_min = p.ai[k_min];
    size_t j_min = MAX(s_j_min, m_j_min);   // столбцы j < j_min нулевые для строки матрицы и строки вектора s
    size_t k = k_max - 1;
    for(;;) // обход элементов строки
    {
        size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
        //if(s[j] != 0)
            E += e.a[k] * s[j] * e.d[j];// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
        if(j <= j_min)
            break;
        k--;
    }
    return E;
}

// вычисление скалярного произведения строки i1 и i2 матрицы (p, e)
// медленно работает
/*double SSCMStrScalarMulSSCMStr(const size_t i1, const size_t i2, const SSCMPortrait &p, const SSCMElements &e)
{
    // слева направо
    size_t k1 = p.ind[i1];          // k1 - индекс элемента строки i1
    size_t k1_max = p.ind[i1 + 1];
    size_t k2 = p.ind[i2];          // k2 - индекс элемента строки i2
    size_t k2_max = p.ind[i2 + 1];
    double E = 0;
    if(k1_max - k1 == 0 ||
       k2_max - k2 == 0)
        return 0;
    size_t j1 = p.ai[k1];  // номер столбца в строке i1
    size_t j2 = p.ai[k2];  // номер столбца в строке i2
    for(;;)
    {
        if(j1 == j2)
        {
            E += e.a[k1] * e.a[k2];
            k1++;
            if(k1 == k1_max)
                return E;
            k2++;
            if(k2 == k2_max)
                return E;
            j1 = p.ai[k1];
            j2 = p.ai[k2];
        }
        else
        {
            if(j1 < j2)
            {
                for(;;)
                {
                    k1++;
                    if(k1 == k1_max)
                        return E;
                    j1 = p.ai[k1];
                    if(j1 >= j2)
                        break;
                }
            }
            else
            {
                for(;;)
                {
                    k2++;
                    if(k2 == k2_max)
                        return E;
                    j2 = p.ai[k2];
                    if(j1 <= j2)
                        break;
                }
            }
        }
    }
    return E;
}*/
//справа налево
/*
size_t k1 = p.ind[i1 + 1];          // k1 - индекс элемента строки i1
size_t k1_min = p.ind[i1];
size_t k2 = p.ind[i2 + 1];          // k2 - индекс элемента строки i2
size_t k2_min = p.ind[i2];
double E = 0;
if(k1 == k1_min ||
   k2 == k2_min)
    return E;
k1--;
k2--;
size_t j1 = p.ai[k1];  // номер столбца в строке i1
size_t j2 = p.ai[k2];  // номер столбца в строке i2

for(;;)
{
    if(j1 == j2)
    {
        E += e.a[k1] * e.a[k2];
        if(k1 == k1_min)
            return E;
        k1--;
        if(k2 == k2_min)
            return E;
        k2--;
        j1 = p.ai[k1];
        j2 = p.ai[k2];
    }
    else
    {
        if(j1 > j2)
        {
            for(;;)
            {
                if(k1 == k1_min)
                    return E;
                k1--;
                j1 = p.ai[k1];
                if(j1 <= j2)
                    break;
            }
        }
        else
        {
            for(;;)
            {
                if(k2 == k2_min)
                    return E;
                k2--;
                j2 = p.ai[k2];
                if(j1 >= j2)
                    break;
            }
        }
    }
}
return E;
*/

SSCMPortraitBulder::SSCMPortraitBulder(const size_t newMatrixSize)
{
    init(newMatrixSize);
}
void SSCMPortraitBulder::init(const size_t newMatrixSize)
{
    m.resize(newMatrixSize);
}
void SSCMPortraitBulder::completePortrait(SlauSolving::SSCMPortrait &p)
{
    size_t matrixSize = m.size();
    size_t elementsNumber = 0;
    for(size_t i = 0; i < matrixSize; i++)
        elementsNumber += m[i].size();
    p.init(elementsNumber, matrixSize);
    size_t count_a = 0;
    for(size_t i = 0; i < matrixSize; i++)
    {
        p.ind[i] = count_a;
        for(std::set<size_t>::iterator it = m[i].begin(); it != m[i].end(); ++it)
        {
            p.ai[count_a] = *it;
            count_a++;
        }
    }
    p.ind[matrixSize] = count_a;
    m.clear();
}

SSCMPreconditioner_base *SSCMPreconditioner_base::gen(const Preconditioning preconditioningType)
{
    switch (preconditioningType)
    {
    case Preconditioning::none:
    {
        return new SSCMPreconditioner_none;
    }break;
    case Preconditioning::diag:
    {
        return new SSCMPreconditioner_diag;
    }break;
    case Preconditioning::SSCM_LLT_0:
    {
        return new SSCMPreconditioner_SSCM_LLT_0;
    }break;
    case Preconditioning::SSCM_LDLT_0:
    {
        return new SSCMPreconditioner_SSCM_LDLT_0;
    }break;
    case Preconditioning::Profile_LLT:
    {
        return new SSCMPreconditioner_Profile_LLT;
    }break;
    case Preconditioning::Profile_LDLT:
    {
        return new SSCMPreconditioner_Profile_LDLT;
    }break;
    case Preconditioning::SSCM_LLT:
    {
        return new SSCMPreconditioner_SSCM_LLT;
    }break;
    case Preconditioning::SSCM_LDLT:
    {
        return new SSCMPreconditioner_SSCM_LDLT;
    }break;
    case Preconditioning::PARDISO:
    {
        return nullptr;
    }break;
    }
}
void SSCMPreconditioner_base::bulid(const SSCM &matrix, const size_t firstStr, double &time)
{
    double time1;
    initPortraitAndAllocateMemory(*matrix.p, time1);
    double time2;
    updatePreconditioner(*matrix.e, firstStr, time2);
    time = time1 + time2;
}
void SSCMPreconditioner_base::InversedMmulVector(Vector &x) const
{
    InversedSmulVector(x);
    InversedQmulVector(x);
}

void SSCMPreconditioner_none::initPortraitAndAllocateMemory(const SSCMPortrait &, double &time)
{
fprintf(stderr, "SSCMPreconditioner_none initPortraitAndAllocateMemory..\n");
    time = 0;
}
void SSCMPreconditioner_none::updatePreconditioner(const SSCMElements &, const size_t, double &time)
{
    time = 0;
}
void SSCMPreconditioner_none::release()
{
}
void SSCMPreconditioner_none::saveBMP(const char *, size_t) const
{
}
void SSCMPreconditioner_none::saveProperties(const char *) const
{
}
size_t SSCMPreconditioner_none::getElementsNumber() const
{
    return 0;
}
void SSCMPreconditioner_none::InversedSmulVector(Vector &) const
{
}
void SSCMPreconditioner_none::InversedQmulVector(Vector &) const
{
}

void SSCMPreconditioner_diag::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_diag initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    d.resize(p0.matrixSize);
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_diag::updatePreconditioner(const SSCMElements &set_e, const size_t, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = d.size();
    for (size_t i = 0; i < matrixSize; i++)
        d[i] = 1./e0.d[i];
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_diag::release()
{
    d.clear();
}
void SSCMPreconditioner_diag::saveBMP(const char *, size_t) const
{
}
void SSCMPreconditioner_diag::saveProperties(const char *fileName) const
{
    FILE *f = fopen(fileName, "w");
    fprintf(f, "Size = %dx\n", (int)d.size());
    fclose(f);
}
size_t SSCMPreconditioner_diag::getElementsNumber() const
{
    return d.size();
}
void SSCMPreconditioner_diag::InversedSmulVector(Vector &x) const
{
    // (d^-1)*x -> x
    const size_t matrixSize = d.size();
    for (size_t i = 0; i < matrixSize; i++)
            x[i] *= d[i];
}
void SSCMPreconditioner_diag::InversedQmulVector(Vector &) const
{

}

void SSCMPreconditioner_SSCM_LLT_0::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_SSCM_LLT_0 initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    m.init(p0.elementsNumber, p0.matrixSize);
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_SSCM_LLT_0::updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    double E;
    for (size_t j = firstStr; j < matrixSize; j++)               // j - номер столбца
    {                                       // i - номер строки
        size_t i;
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && (i = p0.ai[k]) < j; k++)
        {
            size_t k1 = p0.ind[j];
            size_t k2 = p0.ind[i];
            E = 0;
            size_t i1;
            size_t i2;
            for (; k1 < p0.ind[j + 1] && k2 < p0.ind[i + 1] && (i1 = p0.ai[k1]) < i && (i2 = p0.ai[k2]) < i;)
            {
                if (i1 == i2)               // L(j, i1) = L[k1]
                    E -= m.a[k1++] * m.a[k2++]; // L(i, i2) = L[k2]
                else
                    if (i1 < i2)
                        k1++;
                    else
                        k2++;
            }
            m.a[k] = (E + e0.a[k]) / m.d[i];       // L(j, i) = L[k]
        }
        E = 0;
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && p0.ai[k] < j; k++)
            E -= m.a[k] * m.a[k];           // i = ai[k] - номер строки
        m.d[j] = sqrt(E + e0.d[j]);          // A(j, j) = d[j]
                                            // L(j, j) = D[j]
    }
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_SSCM_LLT_0::release()
{
    m.release();
}
void SSCMPreconditioner_SSCM_LLT_0::saveBMP(const char *fileName, size_t size) const
{
    p_ptr0_SSCM->saveBMP(fileName, size);
}
void SSCMPreconditioner_SSCM_LLT_0::saveProperties(const char *fileName) const
{
    p_ptr0_SSCM->saveProperties(fileName);
}
size_t SSCMPreconditioner_SSCM_LLT_0::getElementsNumber() const
{
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    return p0.getElementsNumber();
}
void SSCMPreconditioner_SSCM_LLT_0::InversedSmulVector(Vector &x) const
{
    // (L^-1)*x -> x
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t matrixSize = p.matrixSize;
    double E;
    for (size_t j = 0; j < matrixSize; j++) // j - номер столбца
    {
        E = x[j];
        for (size_t k = p.ind[j]; k < p.ind[j + 1]; k++)
            E -= m.a[k] * x[p.ai[k]]; // L(j, i) = L[k], i = p.ai[k] - номер строки
        x[j] = E / m.d[j];
    }
}
void SSCMPreconditioner_SSCM_LLT_0::InversedQmulVector(Vector &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t j = p.matrixSize - 1;
    for (;;)            // j - номер столбца
    {
        x[j] /= m.d[j];                       // i = ai[k] - номер строки
        for(size_t k = p.ind[j]; k < p.ind[j + 1]; k++)
            x[p.ai[k]] -= x[j] * m.a[k];        // U(i, j) = L[k]
        if(j == 0)
            break;
        j--;
    }
}

void SSCMPreconditioner_SSCM_LDLT_0::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_SSCM_LDLT_0 initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    m.init(p0.elementsNumber, p0.matrixSize);
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_SSCM_LDLT_0::updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    double E;
    for (size_t j = firstStr; j < matrixSize; j++)               // j - номер столбца
    {                                       // i - номер строки
        size_t i;
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && (i = p0.ai[k]) < j; k++)
        {
            size_t k1 = p0.ind[j];
            size_t k2 = p0.ind[i];
            E = 0;
            size_t i1;
            size_t i2;
            for (; k1 < p0.ind[j + 1] && k2 < p0.ind[i + 1] && (i1 = p0.ai[k1]) < i && (i2 = p0.ai[k2]) < i;)
            {                               // L(j, i1) = L[k1]
                if (i1 == i2)               // L(i, i2) = L[k2]
                    E -= m.a[k1++] * m.a[k2++] * m.d[i1];
                else
                    if (i1 < i2)
                        k1++;
                    else
                        k2++;
            }
            m.a[k] = (E + e0.a[k]) / m.d[i];       // L(j, i) = L[k]
        }
        E = 0;
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && (i = p0.ai[k]) < j; k++)
            E -= m.a[k] * m.a[k] * m.d[i];  // L(j, i) = L[k]
                                            // L(i, i) = D[i]
        m.d[j] = E + e0.d[j];                    // A(i, i) = d[i]
    }
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_SSCM_LDLT_0::release()
{
    m.release();
}
void SSCMPreconditioner_SSCM_LDLT_0::saveBMP(const char *fileName, size_t size) const
{
    p_ptr0_SSCM->saveBMP(fileName, size);
}
void SSCMPreconditioner_SSCM_LDLT_0::saveProperties(const char *fileName) const
{
    p_ptr0_SSCM->saveProperties(fileName);
}
size_t SSCMPreconditioner_SSCM_LDLT_0::getElementsNumber() const
{
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    return p0.getElementsNumber();
}
void SSCMPreconditioner_SSCM_LDLT_0::InversedSmulVector(Vector &x) const
{
    // (D^-1*L^-1)*x -> x
    // LD*x_new = x
    // LD(i, j) = L(i, j)*D(j, j)
    // LD(i, i) = D(i, i)
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t matrixSize = p.matrixSize;
    double E;
    for (size_t j = 0; j < matrixSize; j++) // j - номер столбца
    {
        E = x[j];                           // i - номер строки
        for (size_t k = p.ind[j]; k < p.ind[j + 1]; k++)
            E -= m.a[k] * m.d[p.ai[k]] * x[p.ai[k]]; // L(j, i) = L[k], i = p.ai[k] - номер строки
        x[j] = E / m.d[j];
    }
}
void SSCMPreconditioner_SSCM_LDLT_0::InversedQmulVector(Vector &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t j = p.matrixSize - 1;
    for(;;)                    // j - номер столбца
    {
        for(size_t k = p.ind[j]; k < p.ind[j + 1]; k++)  // i = ai[k] - номер строки
            x[p.ai[k]] -= x[j] * m.a[k];                // U(i, j) = L[k]
        if(j == 0)
            break;
        j--;
    }
    /*
    for (j = (int)p.matrixSize - 1; j >= 0; j--)                    // j - номер столбца
        for (k = p.ind[j + 1] - 1; k >= p.ind[j]; k--)  // i = ai[k] - номер строки
            x[p.ai[k]] -= x[j] * m.a[k];                // U(i, j) = L[k]

    */
}

void SSCMPreconditioner_Profile_LLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_Profile_LLT initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    // инициализация портрета
    p.init_LLT(p0);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    // матрицы
    e.init(p.getElementsNumber(), p.getMatrixSize());
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_Profile_LLT::updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // копирование исходной матрицы в профильный формат
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCM3x3StrCopyToProfileStr(i, p, e, p0, e0); // копирование строки
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
int lastPerc = 0;
double lastTime = 0;
    // построение матрицы предобусловливателя
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        const size_t k1_min = p.ind[i];
        const size_t k1_max = p.ind[i + 1]; // строго меньше
        const size_t j1_min = i - (k1_max - k1_min);//p.getFirstCol(i);
        double E_sqr_a = 0;
        size_t k = k1_min;
        size_t j = j1_min; // j - номер стоблца
        for (;k < k1_max;)
        {
            const size_t k2_min = p.ind[j];
            const size_t k2_max = p.ind[j + 1];
            size_t k1 = k1_min;
            size_t k2 = k2_min;
            size_t j2_min = j - (k2_max - k2_min);
            if(j2_min >= j1_min)
            {
                k1 += j2_min - j1_min;
            }
            else
            {
                k2 -= j2_min - j1_min;
            }
            double E = 0;
            size_t k1_max_for_cycle = k1 + (k2_max - k2);// всегда (k2_max - k2) < (k1_max - k1) т.к. начинаем идти с одного столбца
            while(k1 < k1_max_for_cycle)
            {
                E += e.a[k1] * e.a[k2];
                k1++; // L(i, j1) = m.a[k1]
                k2++; // L(j, j2) = m.a[k2]
            }
            e.a[k] = (e.a[k] - E) / e.d[j]; // L(i, j) = m.a[k]
            E_sqr_a += e.a[k] * e.a[k];
            k++;
            j++;
        }
        e.d[i] = sqrt(e.d[i] - E_sqr_a);
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCMPreconditioner_Profile_LLT::release()
{
    p.release();
    e.release();
}
void SSCMPreconditioner_Profile_LLT::saveBMP(const char *fileName, size_t size) const
{
    p.saveBMP(fileName, size);
}
void SSCMPreconditioner_Profile_LLT::saveProperties(const char *fileName) const
{
    p.saveProperties(fileName);
}
size_t SSCMPreconditioner_Profile_LLT::getElementsNumber() const
{
    return p.getElementsNumber();
}
void SSCMPreconditioner_Profile_LLT::InversedSmulVector(Vector &x) const
{
    size_t matrixSize = p.matrixSize;
    // (L^-1)*x -> x
    /*
    for (int i = 0; i < matrixSize; i++)    // i - номер строки
    {
        double E = x[i];
        int j = p.getFirstCol(i);             // j - номер столбца
        for (int k = p.ind[i]; k < p.ind[i + 1] && j < i; k++)
        {
            E -= m.a[k] * x[j];               // L(i, j) = m.a[k]
            j++;
        }
        x[i] = E / m.d[i];
    }
    */
    for (size_t i = 0; i < matrixSize; i++)    // i - номер строки
    {
        double E = 0;
        size_t k = p.ind[i];
        size_t j = p.getFirstCol(i);
        size_t k_max = p.ind[i + 1];
        while(k < k_max)
        {
            E += e.a[k] * x[j];               // L(i, j) = m.a[k]
            k++;
            j++;
        }
        /*
        double E = 0;
        int k_min = p.ind[i];
        int k_max = p.ind[i + 1];
        int j_min = p.getFirstCol(i);
        int j_max = i;
         int l_k = p.ind[i + 1] - p.ind[i];
         int l_j = i - p.getFirstCol(i); // = i - (i - (ind[i + 1] - ind[i])) = l_k
        int k = k_min;
        int j = j_min;
        for (; k < k_max && j < i; )
        {
            E += m.a[k] * x[j];               // L(i, j) = m.a[k]
            k++;
            j++;
        }
        */
        x[i] = (x[i] - E) / e.d[i];
    }
}
void SSCMPreconditioner_Profile_LLT::InversedQmulVector(Vector &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    size_t j = p.matrixSize - 1;
    for (;;)  // j - номер столбца
    {
        x[j] /= e.d[j];
        const size_t k_min = p.ind[j];
        const size_t k_max = p.ind[j + 1];  // строго меньше
        size_t i = j - (k_max - k_min);
        size_t k = k_min;
        while(k < k_max)
        {
            x[i] -= x[j] * e.a[k];         // U(i, j) = L[k]
            i++;    // i < j
            k++;
        }
        if(j == 0)
            break;
        j--;
    }
    /*
    // не оптимально для узкой ленты (int j = i + 1; j < _Высота_(i); j++)
    for (int i = p.matrixSize - 1; i >= 0; i--)  // i - номер строки
    {
        double E = 0;
        for (int j = i + 1; j < p.matrixSize; j++)    // j - номер столбца
        {                                             // L(i,j) = L(j,i) => j - строка, i - столбец (для L)
            int k = p.ind[j + 1] - (j - i);
            if(k >= p.ind[j])
                E += m.a[k]*x[j];                     // m.a[k] = L(j,i)
        }
        x[i] = (x[i] - E) / m.d[i];
    }
    */
}

void SSCMPreconditioner_Profile_LDLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_Profile_LDLT initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    // инициализация портрета
    p.init_LLT(p0);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    // матрицы
    e.init(p.getElementsNumber(), p.getMatrixSize());
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_Profile_LDLT::updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // копирование исходной матрицы в профильный формат
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCM3x3StrCopyToProfileStr(i, p, e, p0, e0); // копирование строки
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
int lastPerc = 0;
double lastTime = 0;
    // построение матрицы предобусловливателя
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        const size_t k1_min = p.ind[i];
        const size_t k1_max = p.ind[i + 1]; // строго меньше
        const size_t j1_min = i - (k1_max - k1_min);//p.getFirstCol(i);
        double E_sqr_a = 0;
        size_t k = k1_min;
        size_t j = j1_min; // j - номер стоблца
        for (;k < k1_max;)
        {
            const size_t k2_min = p.ind[j];
            const size_t k2_max = p.ind[j + 1];
            size_t k1 = k1_min;
            size_t k2 = k2_min;
            size_t j2_min = j - (k2_max - k2_min);
            size_t j_count; // j_count = max(j1_min, j2_min) - номер столбца, одинаковый для перемножаемых строк i и j
            if(j2_min >= j1_min)
            {
                k1 += j2_min - j1_min;
                j_count = j2_min;
            }
            else
            {
                k2 -= j2_min - j1_min;
                j_count = j1_min;
            }
            double E = 0;
            size_t k1_max_for_cycle = k1 + (k2_max - k2);// всегда (k2_max - k2) < (k1_max - k1) т.к. начинаем идти с одного столбца
            while(k1 < k1_max_for_cycle)
            {
                E += e.a[k1] * e.a[k2] * e.d[j_count];
                k1++; // L(i, j_count) = e.a[k1]
                k2++; // L(j, j_count) = e.a[k2]
                j_count++;
            }
            e.a[k] = (e.a[k] - E) / e.d[j]; // L(i, j) = e.a[k]
            E_sqr_a += e.a[k] * e.a[k] * e.d[j];
            k++;
            j++;
        }
        e.d[i] -= E_sqr_a;
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCMPreconditioner_Profile_LDLT::release()
{
    p.release();
    e.release();
}
void SSCMPreconditioner_Profile_LDLT::saveBMP(const char *fileName, size_t size) const
{
    p.saveBMP(fileName, size);
}
void SSCMPreconditioner_Profile_LDLT::saveProperties(const char *fileName) const
{
    p.saveProperties(fileName);
}
size_t SSCMPreconditioner_Profile_LDLT::getElementsNumber() const
{
    return p.getElementsNumber();
}
void SSCMPreconditioner_Profile_LDLT::InversedSmulVector(Vector &x) const
{
    size_t matrixSize = p.matrixSize;
    // (D^-1*L^-1)*x -> x
    // LD*x_new = x
    // LD(i, j) = L(i, j)*D(j, j)
    // LD(i, i) = D(i, i)
    for (size_t i = 0; i < matrixSize; i++)
    {
        double E = 0;
        size_t k = p.ind[i];
        size_t j = p.getFirstCol(i);
        size_t k_max = p.ind[i + 1];
        while(k < k_max)
        {
            E += e.a[k] * e.d[j] * x[j]; // L(i, j) = m.a[k]
            k++;
            j++;
        }
        x[i] = (x[i] - E) / e.d[i];
    }
}
void SSCMPreconditioner_Profile_LDLT::InversedQmulVector(Vector &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    size_t j = p.matrixSize - 1;
    for (;;) // j - номер столбца
    {
        const size_t k_min = p.ind[j];
        const size_t k_max = p.ind[j + 1]; // строго меньше
        size_t i = j - (k_max - k_min);
        size_t k = k_min;
        while(k < k_max)
        {
            x[i] -= x[j] * e.a[k]; // U(i, j) = L[k]
            i++;    // i < j
            k++;
        }
        if(j == 0)
            break;
        j--;
    }
}

void SSCMPreconditioner_SSCM_LLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_SSCM_LLT initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // p0, e0 - портрет и элементы исходной матрицы в строчно-столбцовом формате
    // требуется построить портрет p и элементы m матрицы L
    SSCMPortrait &p = p_SSCM;
    SSCMElements &m = m_SSCM;
    p.init_LLT(p0);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    m.init(p.elementsNumber, matrixSize);
    s.resize(matrixSize);
    for (size_t i = 0; i < matrixSize; i++)
        s[i] = 0;
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_SSCM_LLT::updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    SSCMPortrait &p = p_SSCM;
    SSCMElements &m = m_SSCM;
    // 2) копирование исходной матрицы в L
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCMStrCopyToSSCMStr(i, p, m, p0, e0);    // p0,e -> p,m
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
int lastPerc = 0;
double lastTime = 0;
    // 3) построение разложения - заполнение элементов m
    for (size_t i = firstStr; i < matrixSize; i++)           // i - номер строки
    {
        const size_t k_min = p.ind[i];
        const size_t k_max = p.ind[i + 1];  // строго меньше
        double E_sqr_a = 0;
        size_t k = k_min;
        while(k < k_max)
        {
            size_t j = p.ai[k]; // номер столбца в строке i
            // вычисление скалярного произведения
            size_t s_j_min = p.ai[k_min];
            double E = SSCMStrScalarMulVector(j, p, m, s, s_j_min);
            double newValue = (m.a[k] - E) / m.d[j];
            m.a[k] = newValue; // L(i, j) = m.a[k]
            E_sqr_a += newValue * newValue;
            s[j] = newValue;
            k++;
            /*
            double E = SSCMStrScalarMulSSCMStr(i, j, p, m);
            m.a[k] = (m.a[k] - E) / m.d[j]; // L(i, j) = m.a[k]
            E_sqr_a += m.a[k] * m.a[k];
            k++;
            */
        }
        // очистка всех ненулевых элементов временного вектора
        k = k_min;
        while(k < k_max)
        {
            s[p.ai[k]] = 0;
            k++;
        }
        m.d[i] = sqrt(m.d[i] - E_sqr_a);
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCMPreconditioner_SSCM_LLT::release()
{
    p_SSCM.release();
    m_SSCM.release();
    s.clear();
}
void SSCMPreconditioner_SSCM_LLT::saveBMP(const char *fileName, size_t size) const
{
    p_SSCM.saveBMP(fileName, size);
}
void SSCMPreconditioner_SSCM_LLT::saveProperties(const char *fileName) const
{
    p_SSCM.saveProperties(fileName);
}
size_t SSCMPreconditioner_SSCM_LLT::getElementsNumber() const
{
    return p_SSCM.getElementsNumber();
}
void SSCMPreconditioner_SSCM_LLT::InversedSmulVector(Vector &x) const
{
    // (L^-1)*x -> x
    const SSCMPortrait &p = p_SSCM;
    const SSCMElements &m = m_SSCM;
    const size_t matrixSize = p.matrixSize;
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        const size_t k_min = p.ind[i];
        const size_t k_max = p.ind[i + 1];  // строго меньше
        double E = 0;
        size_t k = k_min;
        while(k < k_max)
        {
            size_t j = p.ai[k]; // номер столбца в строке i
            E += m.a[k] * x[j]; // L(i, j) = m.a[k]
            k++;
        }
        x[i] = (x[i] - E) / m.d[i];
    }
}
void SSCMPreconditioner_SSCM_LLT::InversedQmulVector(Vector &x) const
{
    // (U^-1)*x -> x
    const SSCMPortrait &p = p_SSCM;
    const SSCMElements &m = m_SSCM;
    size_t j = p.matrixSize - 1;
    for (;;)  // j - номер столбца
    {
        x[j] /= m.d[j];
        const size_t k_min = p.ind[j];
        const size_t k_max = p.ind[j + 1];  // строго меньше
        size_t k = k_min;
        while(k < k_max)
        {
            size_t i = p.ai[k]; // номер строки
            x[i] -= x[j] * m.a[k]; // U(i, j) = L[k]
            k++;
        }
        if(j == 0)
            break;
        j--;
    }
}

void SSCMPreconditioner_SSCM_LDLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner_SSCM_LDLT initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // p0, e0 - портрет и элементы исходной матрицы в строчно-столбцовом формате
    // требуется построить портрет p и элементы m матрицы L
    SSCMPortrait &p = p_SSCM;
    SSCMElements &m = m_SSCM;
    p.init_LLT(p0);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    m.init(p.elementsNumber, matrixSize);
    s.resize(matrixSize);
    for (size_t i = 0; i < matrixSize; i++)
        s[i] = 0;
    time = t.getCurrentDuration();
}
void SSCMPreconditioner_SSCM_LDLT::updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    SSCMPortrait &p = p_SSCM;
    SSCMElements &m = m_SSCM;
    // 2) копирование исходной матрицы в L
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCMStrCopyToSSCMStr(i, p, m, p0, e0);    // p0,e -> p,m
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
int lastPerc = 0;
double lastTime = 0;
    // 3) построение разложения - заполнение элементов m
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        const size_t k_min = p.ind[i];
        const size_t k_max = p.ind[i + 1];  // строго меньше
        double E_lld = 0;
        size_t k = k_min;
        while(k < k_max)
        {
            size_t j = p.ai[k]; // номер столбца в строке i
            // вычисление скалярного произведения
            size_t s_j_min = p.ai[k_min];
            double E = SSCMStrScalarMulVectoMulD(j, p, m, s, s_j_min);
            double newValue = (m.a[k] - E) / m.d[j];
            m.a[k] = newValue; // L(i, j) = m.a[k]
            E_lld += newValue * newValue * m.d[j];
            s[j] = newValue;
            k++;
        }
        // очистка всех ненулевых элементов временного вектора
        k = k_min;
        while(k < k_max)
        {
            s[p.ai[k]] = 0;
            k++;
        }
        m.d[i] -= E_lld;
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCMPreconditioner_SSCM_LDLT::release()
{
    p_SSCM.release();
    m_SSCM.release();
    s.clear();
}
void SSCMPreconditioner_SSCM_LDLT::saveBMP(const char *fileName, size_t size) const
{
    p_SSCM.saveBMP(fileName, size);
}
void SSCMPreconditioner_SSCM_LDLT::saveProperties(const char *fileName) const
{
    p_SSCM.saveProperties(fileName);
}
size_t SSCMPreconditioner_SSCM_LDLT::getElementsNumber() const
{
    return p_SSCM.getElementsNumber();
}
void SSCMPreconditioner_SSCM_LDLT::InversedSmulVector(Vector &x) const
{
    // (D^-1*L^-1)*x -> x
    // LD*x_new = x
    // LD(i, j) = L(i, j)*D(j, j)
    // LD(i, i) = D(i, i)
    const SSCMPortrait &p = p_SSCM;
    const SSCMElements &m = m_SSCM;
    size_t matrixSize = p.matrixSize;
    double E;
    for (size_t j = 0; j < matrixSize; j++) // j - номер столбца
    {
        E = x[j]; // i - номер строки
        for (size_t k = p.ind[j]; k < p.ind[j + 1]; k++)
            E -= m.a[k] * m.d[p.ai[k]] * x[p.ai[k]]; // L(j, i) = L[k], i = p.ai[k] - номер строки
        x[j] = E / m.d[j];
    }
}
void SSCMPreconditioner_SSCM_LDLT::InversedQmulVector(Vector &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    const SSCMPortrait &p = p_SSCM;
    const SSCMElements &m = m_SSCM;
    size_t j = p.matrixSize - 1;
    for(;;)                    // j - номер столбца
    {
        for(size_t k = p.ind[j]; k < p.ind[j + 1]; k++)  // i = ai[k] - номер строки
            x[p.ai[k]] -= x[j] * m.a[k];                // U(i, j) = L[k]
        if(j == 0)
            break;
        j--;
    }
}



SSCMSolver_base *SSCMSolver_base::gen(const SolverType solverType)
{
    switch (solverType)
    {
    case SolverType::Direct:
    {
        return new SSCMSolver_direct;
    }break;
    case SolverType::LOS:
    {
        return new SSCMSolver_LOS;
    }break;
    case SolverType::CGM:
    {
        return new SSCMSolver_CGM;
    }break;
    }
}

void SSCMSolver_direct::init(const size_t matrixSize)
{
    r.resize(matrixSize);
    r0.resize(matrixSize);
}
void SSCMSolver_direct::release()
{
    r.clear();
    r0.clear();
}
void SSCMSolver_direct::solve(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, const Vector &x0, const SolverParameters &, Vector &x, double &residual, double &relativeResidual, int &iterations, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    double norma_b = sqrt(VectorScalMul(b, b));
    if(norma_b < 1.e-20)
    {
        iterations = 0;
        relativeResidual = -1;
        residual = 0;
        for(size_t i = 0; i < x.size(); i++)
            x[i] = 0;
        goto ext;
    }
    // прямой решатель
    iterations = 0;
    // невязка первого приближения
    // r0 = A-bx0
    SSCMsolveResidual(matrix, b, x0, r0);
    // решение СЛАУ
    // x = Mm1 * b
    VectorCopy(b, x);  // b -> x
    p->InversedMmulVector(x);
    // невязка решения
    // r = A-bx
    SSCMsolveResidual(matrix, b, x, r);
    // residual = |b-Ax|/|b|
    residual = sqrt(VectorScalMul(r, r)) / norma_b;

    if(VectorScalMul(r0, r0) == 0)
        relativeResidual = 100000;
    else
        relativeResidual = sqrt(VectorScalMul(r, r) / VectorScalMul(r0, r0));
ext:;
    time = t.getCurrentDuration();
}

void SSCMSolver_iterative_base::solve(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, const Vector &x0, const SolverParameters &ssp, Vector &x, double &residual, double &relativeResidual, int &iterations, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    TimeIntervals::timeInterval debug_t;
    debug_t.begin();
    double norma_b = sqrt(VectorScalMul(b, b));

    if(norma_b < 1.e-20)
    {
        iterations = 0;
        relativeResidual = 1;
        residual = 0;
        for(size_t i = 0; i < x.size(); i++)
            x[i] = 0;
        goto ext;
    }
    VectorCopy(x0, x);
    // r0 = A-bx0
    SSCMsolveResidual(matrix, b, x, r0);
    residual = sqrt(VectorScalMul(r0, r0)) / norma_b;
    iterations = 0;
    // итерации
    for(;;)
    {
        //printf("residual = %le\n", residual);
        if (iterations >= ssp.maxIter)
            break;
        if (residual <= ssp.eps)
            break;
        /*
        if (residual <= ssp.eps &&
            residual >= lastResidual)     // дожимаем
            break;*/
        iterations++;
        //if(iterations % 16 == 0)    // периодические обновления процесса
        //    start(portrait, e, b, p, x);
        start(matrix, b, p, x);
        iter(matrix, b, p, x);
        //lastResidual = residual;
        residual = sqrt(VectorScalMul(r, r)) / norma_b;
        //residual = sqrt(r_mul_r) / norma_F;
        if (debug_t.getCurrentDuration() >= 60)   // сообщаем каждую минуту
        {
            //printf("\niter = %6d	residual = %le\n", iterations, residual);
            debug_t.begin();
        }
        //printf("\niter = %6d	residual = %le\n", iterations, residual);
    }
    // r = A-bx
    SSCMsolveResidual(matrix, b, x, r);
    // residual = |b-Ax|/|b|
    residual = sqrt(VectorScalMul(r, r)) / norma_b;
    // relativeResidual = |b-Ax|/|b-Ax0| = |r|/|r0|
    if(VectorScalMul(r0, r0) == 0)
        relativeResidual = 100000;
    else
        relativeResidual = sqrt(VectorScalMul(r, r) / VectorScalMul(r0, r0));
ext:;
    time = t.getCurrentDuration();
}

void SSCMSolver_LOS::init(const size_t matrixSize)
{
    r.resize(matrixSize);
    r0.resize(matrixSize);
    ss.resize(matrixSize);
    z.resize(matrixSize);
    w.resize(matrixSize);
    aa.resize(matrixSize);
    pp.resize(matrixSize);
}
void SSCMSolver_LOS::release()
{
    r.clear();
    r0.clear();
    ss.clear();
    z.clear();
    w.clear();
    aa.clear();
    pp.clear();
}
void SSCMSolver_LOS::start(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, Vector &x)
{
    // r0 = (b - A*x0)
    SSCMmulVector(matrix, x, r);        //A*x -> r
    Vector1PlusCmulVector2(b, -1, r, r);//b - r -> r
    // p0 = s0 = Mm1*r0
    VectorCopy(r, pp);  // r -> p
    p->InversedMmulVector(pp);
    VectorCopy(pp, ss); // p -> ss
    // z0 = aa0 = A*p0
    SSCMmulVector(matrix, pp, z);   // A*p -> z
    VectorCopy(z, aa);  // z -> aa
    // w0 = Mm1*z0
    VectorCopy(z, w);   // z -> w
    p->InversedMmulVector(w);
}
void SSCMSolver_LOS::iter(const SSCM &matrix, const Vector &, const SSCMPreconditioner_base *p, Vector &x)
{
    double ca = VectorScalMul(w, r) / VectorScalMul(w, z);
    Vector1PlusCmulVector2(x, ca, pp, x);   // x + ca * p -> x
    Vector1PlusCmulVector2(r, -ca, z, r);   // r - ca * z -> r
    Vector1PlusCmulVector2(ss, -ca, w, ss); // ss - ca * w -> ss
    SSCMmulVector(matrix, ss, aa);          // A*ss -> aa
    double cb = -(VectorScalMul(w, aa) / VectorScalMul(w, z));
    Vector1PlusCmulVector2(ss, cb, pp, pp); // s + cb * p -> p
    Vector1PlusCmulVector2(aa, cb, z, z);   // aa + cb * z -> z
    // w = Mm1*z
    VectorCopy(z, w);   // z -> w
    p->InversedMmulVector(w);
}

void SSCMSolver_CGM::init(const size_t matrixSize)
{
    r.resize(matrixSize);
    r0.resize(matrixSize);
    ss.resize(matrixSize);
    z.resize(matrixSize);
    w.resize(matrixSize);
    aa.resize(matrixSize);
    pp.resize(matrixSize);
}
void SSCMSolver_CGM::release()
{
    r.clear();
    r0.clear();
    ss.clear();
    z.clear();
    w.clear();
    aa.clear();
    pp.clear();
}
void SSCMSolver_CGM::start(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, Vector &x)
{
    // r0 = (b - A*x0)
    SSCMmulVector(matrix, x, r);        //A*x -> r
    Vector1PlusCmulVector2(b, -1, r, r);//b - r -> r
    // z0 = Mm1*r0
    VectorCopy(r, z);   // r -> z
    p->InversedMmulVector(z);
}
void SSCMSolver_CGM::iter(const SSCM &matrix, const Vector &, const SSCMPreconditioner_base *p, Vector &x)
{
    // aa = Mm1*r
    VectorCopy(r, aa);  // r -> aa
    p->InversedMmulVector(aa);
    // A*z -> pp
    SSCMmulVector(matrix, z, pp);
    double ca1 = VectorScalMul(aa, r);
    double ca = ca1 / VectorScalMul(pp, z);
    // x + ca * z -> x
    Vector1PlusCmulVector2(x, ca, z, x);
    // r - ca * (A*z) -> r
    Vector1PlusCmulVector2(r, -ca, pp, r);
    // ss = Mm1*r
    VectorCopy(r, ss);  // r -> ss
    p->InversedMmulVector(ss);
    double cb = VectorScalMul(ss, r) / ca1;
    // ss + cb*z -> z
    Vector1PlusCmulVector2(ss, cb, z, z);
}

SSCMBulder::SSCMBulder(const size_t matrixSize)
{
    str.resize(matrixSize);
}
SSCMBulder::~SSCMBulder()
{
    release();
}
void SSCMBulder::release()
{
    str.clear();
}
void SSCMBulder::setMatrixSize(const size_t matrixSize)
{
    str.resize(matrixSize);
}
void SSCMBulder::start()
{
    srand(0);   //##для отладки
    for(size_t i = 0; i < str.size(); i++)
    {
        str[i].resize(0);
    }
}
void SSCMBulder::complete(SSCM &matrix)
{
    SSCMPortrait &p = *matrix.p;
    SSCMElements &e = *matrix.e;
    size_t matrixSize = str.size();
    size_t elementsNumber = 0;
    for(size_t i = 0; i < matrixSize; i++)
    {
        // сортировка строки
        if(str[i].size() != 0)
        {
            qsort(&str[i][0], str[i].size(), sizeof(SlauElement), elementsCmp);
            // подсчет различных внедиагональных ячеек матрицы
            size_t jLast = -1;
            for (size_t j = 0; j < str[i].size(); j++)
            {
                if(str[i][j].j == i) break;
                if(str[i][j].j != jLast)
                {
                    elementsNumber++;
                    jLast = str[i][j].j;
                }
            }
        }
    }
    p.init(elementsNumber, matrixSize);
    e.init(elementsNumber, matrixSize);
    // очистка диагонали
    for (size_t i = 0; i < matrixSize; i++)
        e.d[i] = 0;
    // теперь можно представить матрицу в разреженном формате
    // a[ind[i]] - первый элемент строки i. (k = ind[i]...(ind[i+1]-1) для строки i)
    // ai[k] - номер столбца элемента a[k]
    // элементы нижнего треугольника (доступ по индексу k)
    // заполнение ai[] и a[]
    size_t count_a = 0;    // индекс в массивах ai[] и a[]
    p.ind[0] = 0;
    size_t i = 0;
    size_t str_ind = 0;
    for (;;)
    {
        // ищем первый элемент перед началом суммирования
        for(;;)
        {
            if(str_ind >= str[i].size())
            {
                i++;
                str_ind = 0;
                p.ind[i] = count_a;
                if(i >= matrixSize) break;
            }
            else
                break;
        }
        if(i >= matrixSize) break;
        size_t key_i = i;
        size_t key_j = str[i][str_ind].j;    // координаты очередной ячейки
        double E = 0;                     // суммарное значение в этой ячейке
        for (;;)
        {
            E += str[i][str_ind].a;
            str_ind++;
            // суммируем элементы, принадлежащие одной ячейке
            if(str_ind >= str[i].size() || str[i][str_ind].j != key_j) break;
        }
        if(key_i != key_j)
        {
            e.a[count_a] = E;
            p.ai[count_a] = key_j;
            count_a++;
        }
        else
            e.d[key_i] = E;
    }
}
void SSCMBulder::completeWithPortrait(SSCM &matrix)
{
    SSCMPortrait &p = *matrix.p;
    SSCMElements &e = *matrix.e;
    size_t matrixSize = str.size();
    for(size_t i = 0; i < matrixSize; i++)
    {
        // сортировка строки
        if(str[i].size() != 0)
            qsort(&str[i][0], str[i].size(), sizeof(SlauElement), elementsCmp);
    }
    e.init(p.elementsNumber, p.matrixSize);
    // очистка диагонали
    for (size_t i = 0; i < p.matrixSize; i++)
        e.d[i] = 0;
    size_t count_a = 0;    // индекс в массивах ai[] и a[]
    size_t i = 0;
    size_t str_ind = 0;
    for (;;)
    {
        // ищем первый элемент перед началом суммирования
        for(;;)
        {
            if(str_ind >= str[i].size())
            {
                i++;
                str_ind = 0;
                if(i >= matrixSize) break;
            }
            else
                break;
        }
        if(i >= matrixSize) break;
        size_t key_i = i;
        size_t key_j = str[i][str_ind].j;    // координаты очередной ячейки
        double E = 0;                     // суммарное значение в этой ячейке
        for (;;)
        {
            E += str[i][str_ind].a;
            str_ind++;
            // суммируем элементы, принадлежащие одной ячейке
            if(str_ind >= str[i].size() || str[i][str_ind].j != key_j) break;
        }
        if(key_i != key_j)
        {
            e.a[count_a] = E;
            count_a++;
        }
        else
            e.d[key_i] = E;
    }
}
void SSCMBulder::fixReservedMemory()
{
    // временный массив строк
    std::vector<std::vector<SlauElement>> t(str.size());
    for(size_t i = 0; i < str.size(); i++)
    {
        t[i].resize(str[i].size());
        for(size_t j = 0; j < str[i].size(); j++)
            t[i][j] = str[i][j];
    }
    str.clear();
    str.resize(t.size());
    for(size_t i = 0; i < t.size(); i++)
    {
        str[i].resize(t[i].size());
        str[i].reserve(t[i].size());
        for(size_t j = 0; j < t[i].size(); j++)
            str[i][j] = t[i][j];
    }
    t.clear();
    /*
    // резервируем память для следующего сбора матрицы
    for(size_t i = 0; i < str.size(); i++)
        str[i].reserve(str[i].size());
        */
}
int SSCMBulder::elementsCmp(const void* x1, const void* x2)
{
    SlauElement &px1 = *(SlauElement *)x1;
    SlauElement &px2 = *(SlauElement *)x2;
    if (px1.j < px2.j) return -1;
    if (px1.j > px2.j) return 1;
    //if (px1.a > px2.a) return -1;
    //if (px1.a < px2.a) return 1;  // значения тоже сортируются, по убыванию
    double a1 = fabs(px1.a);
    double a2 = fabs(px2.a);
    if (a1 > a2) return -1;
    if (a1 < a2) return 1;  // значения тоже сортируются, по убыванию модуля
    return 0;
}


}   // namespace SlauSolving
