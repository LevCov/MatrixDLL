//
//  MatrixDLL.cpp
//  MatrixDLL
//
//  Created by Лев Шеховцов on 11.06.2022.
//

#include <iostream>
#include "/Users/levsehovcov/Desktop/tst/MatrixDLL/MatrixDLL/MatrixDLL.hpp"






matrix::matrix(std::vector<std::vector<double>> vecVec): matr(vecVec){ }
identity::identity(long N)
{
    matr = {{}};
    if( N != 0)
    {
        matr.resize(N);
        for(long i = 0; i < N; ++i)
        {
            matr[i].resize(N);
        }
        for(long i = 0; i < N; ++i)
        {
            for(long j = 0; j < N; ++j)
            {
                matr[i][j] = 0;
            }
            matr[i][i] = 1;
        }
    }
}

diagonal::diagonal(long N, double* mass)
{
    matr = {{}};
    if( N != 0)
    {
        matr.resize(N);
        for(long i = 0; i < N; ++i)
        {
            matr[i].resize(N);
        }
        for(long i = 0; i < N; ++i)
        {
            for(long j = 0; j < N; ++j)
            {
                matr[i][j] = 0;
            }
            matr[i][i] = mass[i];
        }
    }
}

upper::upper(std::vector<std::vector<double>> vecVec)
{
    matr = {{}};
    if(vecVec.size() != vecVec[0].size()) throw Exception1();
    if( !vecVec.empty())
    {
        const long N = vecVec.size();
        matr.resize(N);
        for(long i = 0; i < N; ++i)
        {
            matr[i].resize(N);
        }
        for(long i = 0; i < N; ++i)
        {
            for(long j = i; j < N; ++j)
            {
                matr[i][j] = vecVec[i][j];
            }
        }
    }
}

lower::lower(std::vector<std::vector<double>> vecVec)
{
    matr = {{}};
    if(vecVec.size() != vecVec[0].size()) throw Exception1();
    if( !vecVec.empty())
    {
        const long N = vecVec.size();
        matr.resize(N);
        for(long i = 0; i < N; ++i)
        {
            matr[i].resize(N);
        }
        for(long i = 0; i < N; ++i)
        {
            for(long j = 0; j <= i; ++j)
            {
                matr[i][j] = vecVec[i][j];
            }
        }
    }
}

symmetric::symmetric(std::vector<std::vector<double>> vecVec)
{
    matr = {{}};
    if(vecVec.size() != vecVec[0].size()) throw Exception1();
    if( !vecVec.empty())
    {
        const long N = vecVec.size();
        matr.resize(N);
        for(long i = 0; i < N; ++i)
        {
            matr[i].resize(N);
        }
        for(long i = 0; i < N; ++i)
        {
            for(long j = i; j < N; ++j)
            {
                matr[i][j] = vecVec[i][j];
                matr[j][i] = vecVec[i][j];;
                
            }
        }
    }
}

matrix operator +(const matrix& M1, const matrix& M2)
{
    if(M1.matr.size() != M2.matr.size() || M1.matr[0].size() != M2.matr[0].size()) throw Exception1();
    
    long N = M1.matr.size();
    long M = M1.matr[0].size();
    matrix Tmp;
    
    Tmp.matr.resize(N);
    for(long i = 0; i < N; ++i) Tmp.matr[i].resize(M);
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        Tmp.matr[i][j] = M1.matr[i][j] + M2.matr[i][j];
    }
    return Tmp;
}

matrix operator -(const matrix& M1, const matrix& M2)
{
    if(M1.matr.size() != M2.matr.size() || M1.matr[0].size() != M2.matr[0].size()) throw Exception1();
    
    long N = M1.matr.size();
    long M = M1.matr[0].size();
    matrix Tmp;
    
    Tmp.matr.resize(N);
    for(long i = 0; i < N; ++i) Tmp.matr[i].resize(M);
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        Tmp.matr[i][j] = M1.matr[i][j] - M2.matr[i][j];
    }
    return Tmp;
}

matrix operator *(const matrix& M1, const matrix& M2)
{
    if(M1.matr[0].size() != M2.matr.size()) throw Exception1();
    matrix Tmp; // (n, m)x(m,k)
    //(n, 1)x(1, n)
    long M = M1.matr[0].size();
    long N = M1.matr.size();
    long K = M2.matr[0].size();
    Tmp.matr.resize(N);
    
    for(long k = 0; k < N; ++k)
    {
        Tmp.matr[k].resize(K);
    }
    
    for(long n = 0; n < N; ++n)
    {
        for(long k = 0; k < K; ++k)
        {
            for(long m = 0; m < M; ++m)
            {
                if(m == 0) Tmp.matr[n][k] = M1.matr[n][m] * M2.matr[m][k];
                else Tmp.matr[n][k] += M1.matr[n][m] * M2.matr[m][k];
            }
        }
    }
    return Tmp;
}

matrix matrix::Adam(const matrix& M2)
{
    if(this->matr[0].size() != M2.matr[0].size() || this->matr.size() != M2.matr.size()) throw Exception1();
    matrix Tmp;
    Tmp.matr.resize(this->matr.size());
    long M = matr[0].size();
    for(long i = 0; i < this->matr.size(); ++i)
    {
        Tmp.matr[i].resize(M);
    }
    
    for(long i = 0; i < this->matr.size(); ++i)
    {
        for(long j = 0; j < M; ++j)
        {
            Tmp.matr[i][j] = this->matr[i][j]*M2.matr[i][j];
        }
    }
    
    return Tmp;
}

matrix identity::Adam(const matrix& M2)
{
    if(this->matr[0].size() != M2.matr[0].size() || this->matr.size() != M2.matr.size()) throw Exception1();
    matrix Tmp(this->matr);
    for(long i = 0; i < this->matr.size(); ++i )
    {
        Tmp.matr[i][i] = M2.matr[i][i];
    }
    
    return Tmp;
}


matrix matrix::mul(double Num)
{
    matrix Tmp(this->matr);
    for(int n = 0; n < this->matr.size(); ++n)
    {
        for(int m = 0; m < this->matr[0].size(); ++m)
        {
            Tmp.matr[n][m] = Tmp.matr[n][m] * Num;
        }
    }
    return Tmp;
}
//2 лаба ----------------------------------------

double matrix::det()
{
    if(this->matr.size() != this->matr[0].size()) throw Exception1();
    matrix Tmp(this->matr);
    bool fl = false;
    long index[this->matr[0].size()] ;
    
    for(long j = 0; j < this->matr[0].size(); ++j)
    {
        fl = false;
        for(long i = 0; i < this->matr.size(); ++i)
        {
            if(this->matr[i][j] != 0)
            {
                fl =  true;
                index[j] = i;
            }
        }
        if(fl == false)
        {
            break;
        }
    }
    //for(int i = 0; i < this->matr[0].size(); ++i) std::cout<<"k= "<<i<<"index[k] = "<<index[i]<<std::endl;
    if(fl == false) return 0.0;
    
    for(long k = 0; k < this->matr[0].size(); ++k)
    {
        if(k == index[k]) continue;
        for(long j = 0; j< this->matr[0].size(); ++j)
        {
            Tmp.matr[k][j] += Tmp.matr[index[k]][j];
        }
    }
    
    double tmp;
    for (int k = 0; k < this->matr.size() - 1; k++) {
        for (int i = k + 1; i < this->matr.size(); i++) {
            tmp = -Tmp.matr[i][k] / Tmp.matr[k][k];
            for (int j = 0; j < this->matr.size(); j++) {
                Tmp.matr[i][j] += Tmp.matr[k][j] * tmp;
            }
        }
        }
        tmp = 1;
    //for(int i = 0; i < this->matr[0].size(); ++i) std::cout<<"k= "<<i<<"index[k] = "<<index[i]<<std::endl;
    for(int i = 0; i < this->matr.size(); ++i)
    {
        tmp *= Tmp.matr[i][i];
    }

    return tmp;
}

double diagonal::det()
{
    double tmp=1;
    for(int i = 0; i < this->matr.size(); ++i)
    {
        tmp *= this->matr[i][i];
    }
    return tmp;
}

double upper::det()
{
    double tmp=1;
    for(int i = 0; i < this->matr.size(); ++i)
    {
        tmp *= this->matr[i][i];
    }
    return tmp;
}

double lower::det()
{
    double tmp=1;
    for(int i = 0; i < this->matr.size(); ++i)
    {
        tmp *= this->matr[i][i];
    }
    return tmp;
}

double matrix::scalar(const matrix& M)
{
    double tmp = 0;
    if(this->matr.size() == 1 && this->matr[0].size() == M.matr.size() && M.matr[0].size() == 1)
    {
        for(int i = 0; i < this->matr[0].size(); ++i)
        {
            tmp += this->matr[0][i] * M.matr[i][0];
        }
    }
    else throw Exception2();
    return tmp;
}

double matrix::Fronebius()
{
    double tmp = 0;
    
    for(int i = 0; i < this->matr.size(); ++i)
    {
        for(int j = 0; j < this->matr.size(); ++j)
        {
            tmp += this->matr[i][j] * this->matr[i][j];
        }
    }
    return pow(tmp, 1./2);
}

double diagonal::Fronebius()
{
    double tmp = 0;
    
    for(int i = 0; i < this->matr.size(); ++i)
    {
        tmp += this->matr[i][i] * this->matr[i][i];
    }
    return pow(tmp, 1./2);
}

double matrix::norm()
{
    double tmp = 0;
    if(this-> matr.size() == 1)
    {
        for(int i = 0; i < this-> matr[0].size(); ++i) tmp += this-> matr[0][i] * this-> matr[0][i];
        return pow(tmp, 1./2);
    }
    else if(this-> matr[0].size() == 1)
    {
        for(int i = 0; i < this-> matr.size(); ++i) tmp += this-> matr[i][0] * this-> matr[i][0];
        return pow(tmp, 1./2);
    }
    else throw Exception2();
}

double matrix::maxnorm()
{
    double tmp = 0;
    if(this-> matr.size() == 1)
    {
        tmp = this-> matr[0][0];
        for(int i = 1; i < this-> matr[0].size(); ++i) if(tmp < this-> matr[0][i]) tmp = this-> matr[0][i];
        return tmp;
    }
    else if(this-> matr[0].size() == 1)
    {
        tmp = this-> matr[0][0];
        for(int i = 0; i < this-> matr.size(); ++i) if(tmp < this-> matr[i][0]) tmp = this-> matr[i][0];
        return tmp;
    }
    else throw Exception2();
}

// 3 лаба ----------------------------------------------------

matrix matrix::transpose()
{
    matrix Tmp;
    long N = this->matr.size();
    long M = this->matr[0].size();
    Tmp.matr.resize(M);
    for(int i = 0; i < M; ++i)
    {
        Tmp.matr[i].resize(N);
    }
    
    for(int n = 0; n < N; ++n)
    {
        for(int m = 0; m < M; ++m)
        {
            Tmp.matr[m][n] = this->matr[n][m];
        }
    }
    return Tmp;
}
double matrix::angle(const matrix& M)
{
    /*if(!((this-> matr.size() == 1 || this-> matr[0].size() == 1) && (M.matr.size() == 1 || M.matr[0].size() == 1)))
    {
        throw Exception2();
    }*/
    
    double tmp = 0;
    double cosangle = 0;
    //считаем норму M
    if(M.matr.size() == 1)
    {
        for(int i = 0; i < M.matr[0].size(); ++i) tmp += M.matr[0][i] * M.matr[0][i];
        tmp = pow(tmp, 1./2);
    }
    else
    {
        for(int i = 0; i < M.matr.size(); ++i) tmp += M.matr[i][0] * M.matr[i][0];
        tmp = pow(tmp, 1./2);
    }
    // посчитали
    
    if(tmp == 0 || (*this).norm() == 0)
    {
        throw Exception3();
    }
    
    if(this->matr[0].size() == M.matr.size() && M.matr[0].size() == 1)
    {
        //std::cout<<1<<std::endl;
        cosangle = (*this).scalar(M)/ tmp / (*this).norm();
    }
    else if(this->matr[0].size() == M.matr[0].size() && M.matr.size() == 1)
    {
        matrix Vec(M.matr);
        Vec = Vec.transpose();
        //std::cout<<2<<std::endl;
        cosangle = (*this).scalar(Vec)/ tmp / (*this).norm();
    }
    else if(this->matr[0].size() == M.matr[0].size() && M.matr[0].size() == 1 )
    {
        //std::cout<<3<<std::endl;
        cosangle = ((*this).transpose()).scalar(M) / tmp / (*this).norm();
    }
    else if(this->matr.size() == M.matr[0].size() && M.matr.size() == 1 )
    {
        matrix Vec(M.matr);
        //std::cout<<4<<std::endl;
        cosangle = Vec.scalar((*this))/ tmp / (*this).norm();
    }
    return acos(cosangle);
}

long matrix::rank()
{
    matrix Tmp(this->matr);
    if(this->matr.size() < this->matr[0].size()) Tmp = Tmp.transpose();
    
    bool fl = false;
    long N = Tmp.matr.size(); //строки
    long M = Tmp.matr[0].size(); //столбцы
    long index[M] ;
    
    for(long j = 0; j < M; ++j)
    {
        fl = false;
        for(long i = 0; i < N; ++i)
        {
            if(this->matr[i][j] != 0)
            {
                fl =  true;
                index[j] = i;
            }
        }
        if(fl == false)
        {
            break;
        }
    }
    if(fl == false) return 0;
    //return Tmp;
    
    for(long k = 0; k < M; ++k)
    {
        if(k == index[k]) continue;
        for(long j = 0; j< M; ++j)
        {
            Tmp.matr[k][j] += Tmp.matr[index[k]][j];
        }
    }
    //return Tmp;
    double tmp;
    for (long k = 0; k < M; k++)
    {
        for (long i = k + 1; i < N; i++)
        {
            tmp = -Tmp.matr[i][k] / Tmp.matr[k][k];
            for (long j = 0; j < M; j++)
            {
                Tmp.matr[i][j] += Tmp.matr[k][j] * tmp;
            }
        }
        }
        long rank = M;
        for(long i = 0; i < M; ++i)
        {
            if(Tmp.matr[i][i] == 0) --rank;
        }
        return rank;
}


std::vector<std::vector<double>> veccopy(std::vector<std::vector<double>> dest, std::vector<std::vector<double>> vec, long n, long m)
{
    int ki = 0;
    for(int i = 0; i < vec.size(); ++i)
    {
        if(i != n)
        {
            for (int j = 0, kj = 0; j < vec.size(); j++)
            {
                    if (j != m)
                    {
                   dest[ki][kj] = vec[i][j];
                    kj++;
                        }
                    }
                ki++;
        }
        
    }
    //std::cout<<dest[0]
    return dest;
}

matrix matrix::inverse()
{
    if(this->matr.size() != this->matr[0].size()) throw Exception4();
    if((*this).rank() != this->matr[0].size()) throw Exception5();
    
    std::vector<std::vector<double>> tmp;
    long N = this->matr[0].size();
    tmp.resize(N - 1);
    for(int i = 0; i < N - 1; ++i)
    {
        tmp[i].resize(N - 1);
    }
    matrix Tmp(tmp), Inv(this->matr);
    //tmp = veccopy(tmp, this->matr, 0, 0);
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            tmp = veccopy(tmp, this->matr, i, j);
            Tmp.matr  = tmp;
            Inv.matr[i][j] = pow(-1, i + j) * Tmp.det() / (*this).det();
            
        }
    }
    return Inv.transpose();
}

// 4 лаба ------------------------------------------------------------------
std::ostream& operator <<(std::ostream& out, const matrix& M)
{
    for (int i = 0; i < M.matr.size(); i++)
    {
        for (int j = 0; j < M.matr[0].size(); j++)
        {
            out << std::right << std::setw(12) << std::setfill(' ') << M.matr[i][j];
        }
        out << std::endl;
    }
    return out;
}


std::istream& operator>>(std::istream& in, matrix& M)
{
    if (in.fail()) throw Exception6();
    if (in.peek() == EOF) throw Exception7();


    std::streampos pos = 0;
    in.seekg(0, in.end);
    pos = in.tellg();
    in.seekg(0, in.beg);

    std::string line, check;
    std::stringstream instr;
    long rows = 1, cols = 0, maxcols = 0;
    double x;

    getline(in, line);
    while (line.find(',') != std::string::npos)
    {
        line.replace(line.find(','), 1, ".");
    }

    instr << line;
    while (!instr.eof())    // цикл для подсчёта столбцов в первой строке
    {
        instr >> x;
        ++cols;
    }
    maxcols = cols;

    while (1)    // цикл для подсчёта строк в матрице
    {
        if (in.tellg() == pos) break;

        getline(in, line);
        ++rows;
    }

    in.seekg(0, in.beg);
    
    M.matr.resize(rows);
    for(long i = 0; i < rows; ++i) M.matr[i].resize(cols);


    getline(in, line);
    while (line.find(',') != std::string::npos)    // цикл для замены всех "," на "."
    {
        line.replace(line.find(','), 1, ".");
    }

    instr.str(std::string());
    instr.clear();
    
    for (int i = 0; i < rows; ++i)
    {
        cols = 0;
        instr << line;
        while (!instr.eof())
        {
            instr >> check;
            try
            {
                x = stod(check);
            }
            catch (std::invalid_argument)
            {
                throw InputError(std::string("The text file contains a value different from numbers"));
            }
            catch (std::out_of_range)
            {
                throw InputError(std::string("The text file contains a value too great to be read"));
            }
            catch (...)
            {
                throw;
            }

            // проверка того, что мы не выйдем за пределы выделенной памяти
            if ( (i * maxcols + cols) >= (rows * maxcols) ) break;

            M.matr[i][cols] = x;
            ++cols;
        }

        if (cols != maxcols) throw SizeError(std::string("Matrix in text file has different dimensions from those that are given"));

        getline(in, line);
        while (line.find(',') != std::string::npos)    // цикл для замены всех "," на "."
        {
            line.replace(line.find(','), 1, ".");
        }

        instr.str(std::string());
        instr.clear();
    }

    instr.str(std::string());

    return in;
}
//-----------------------
void matrix::readBin(std::istream& in)
{
    if (in.fail()) throw Exception6();
    if (in.peek() == EOF) throw Exception7();

    long sizeMatrix[2];
    double x = 0;
    in.read((char*)&sizeMatrix, sizeof(long) * 2);
    long rows = sizeMatrix[0];
    long cols = sizeMatrix[1];
    this->matr.resize(rows);
    
    for(int i = 0; i < rows; ++i) this->matr[i].resize(cols);
    
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (in.peek() == EOF) throw SizeError(std::string("Matrix in binary file has different dimensions from those that are given"));

            in.read((char *)&x, sizeof(double));
            (*this).matr[i][j] = x;
        }
    }
    return;
}

// метод для записи матрицы в бинарный файл
void matrix::writeBin(std::ostream& out)
{
    unsigned long sizeMatrix[2] = { (*this).matr.size(), (*this).matr[0].size() };
    double x;
    out.write((char*)&sizeMatrix, sizeof(long) * 2);
    for (int i = 0; i < this->matr.size(); ++i)
    {
        for (int j = 0; j < this->matr[0].size(); ++j)
        {
            x = (*this).matr[i][j];
            out.write((char*)&x, sizeof(double));
        }
    }
    return;
}
// Лабораторная работа №5 ||||||||||||||||||||||||||||||||||||||||||||||||||||||

matrix PCA::center()
{
    long N = this->mat.matr[0].size();//столбцы
    long M = this->mat.matr.size();
    
    matrix Tmp = this->mat;
    
    double mid[N];
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        {
            mid[i] += this->mat.matr[j][i];
        }
        mid[i] = mid[i]/M;
    }
    
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        {
            Tmp.matr[j][i] -= mid[i];
        }
    }
    
    return Tmp;
}

matrix PCA::scaling()
{
    long N = this->mat.matr[0].size(); // столбцы
    long M = this->mat.matr.size(); // строки
    
    std::vector<double> vec = {0};
    vec.resize(N);
    
    matrix Tmp  = this->mat;
    
    double mid[N];
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        {
            mid[i] += this->mat.matr[j][i];
        }
        mid[i] = mid[i]/M;
    }
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        {
            vec[i] += (this->mat.matr[j][i] - mid[i]) * (this->mat.matr[j][i] - mid[i]) / (M - 1);
        }
        vec[i] = std::sqrt(vec[i]);
    }
    
    for(long i = 0; i < N; ++i)
    {
        for(long j = 0; j < M; ++j)
        {
            if(vec[i] != 0)
                Tmp.matr[j][i] = (Tmp.matr[j][i] - mid[i])/vec[i];
        }
    }
    return Tmp;
}


std::tuple<matrix, matrix, matrix> PCA::NIPALS(long PC)
{
    
    double eps = 1e-8;
    PCA D((*this).scaling());
    
    matrix E(D.mat.matr);
    matrix t, t_old;
    matrix p, d;
    std::vector<std::vector<double>> Pmatr;
    //Pmatr.resize(PC);
    std::vector<std::vector<double>> Tmatr;
    //Tmatr.resize(PC);
    
    for(long h = 0; h < PC; ++h)
    {
        t = {(E.transpose()).matr[h]};
        t = t.transpose();
        do
        {
            p = ((t.transpose() * E).mul(1 / (t.transpose() * t).matr[0][0])).transpose();
            p = p.mul(1 / p.norm());
            t_old = t;
            t = (E * p).mul(1 / (p.transpose() * p).matr[0][0]);
            d = t_old - t;
        }
        while(d.norm() > eps);
        E = E - t * (p.transpose());
        //std::cout<< t<< std::endl;
        //std::cout<< p<< std::endl;
        Pmatr.push_back(p.transpose().matr[0]);
        Tmatr.push_back(t.transpose().matr[0]);
    }
    matrix P(Pmatr);
    P = P.transpose();
    matrix T(Tmatr);
    T = T.transpose();
    return std::make_tuple(T, P, E);
}

matrix PCA::leverage(long PC)// размах
{
    //if(PC < 1 || PC > this->mat.rank()) PC = 1;
    std::tuple<matrix, matrix, matrix> NIPALS_RESULT = this->NIPALS(PC);
    
    matrix T = std::get<0>(NIPALS_RESULT);
    matrix t;
    std::vector<double> hh;
    matrix tmp;

    for(long h = 0; h < T.matr.size(); ++h)
    {
        t = {T.matr[h]};
        t = t.transpose();
        hh.push_back((t.transpose() * ((T.transpose() * T).inverse()) * t).matr[0][0]);
    }
    
    matrix h({hh});
    h = h.transpose();
    return h;
    
}

matrix PCA::deviation(long PC){ //отклонение
    if(PC < 1 || PC > this->mat.rank()) PC = 1;
    std::tuple<matrix, matrix, matrix> NIPALS_RESULT = this->NIPALS(PC);
    
    matrix E = std::get<2>(NIPALS_RESULT);
    
    std::vector<double> dev;
    
    matrix tmp;
    
    for(long h = 0; h < E.matr.size(); ++h)
    {
        tmp = {E.matr[h]};
        dev.push_back(tmp.norm() * tmp.norm());
    }
    tmp = {dev};
    tmp = tmp.transpose();
    return tmp;
}

std::pair<double, double> PCA::dispersion(long PC)
{
    
    matrix v = this->deviation(PC);
    
    double v0 = 0;
    long I = v.matr.size();
    long J = this->mat.matr[0].size();
    for(int h = 0; h < I; ++h)
    {
        v0 += v.matr[h][0] / I;
    }
    
    double TRV = v0 / J;
    double ERV = 0;
    
    matrix X = this->scaling();
    
    for(int i = 0; i < I; ++i)
    {
        for(int j = 0; j < J; ++j)
            ERV += X.matr[i][j] * X.matr[i][j];
    }
    
    ERV = 1 - I * v0 / ERV;
    
    return std::make_pair(ERV, TRV);
}


