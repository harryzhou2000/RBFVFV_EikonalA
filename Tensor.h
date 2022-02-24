#ifndef _TENSOR_H
#define _TENSOR_H
#include "TypeDefine.h"
#include <list>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
namespace ScalarCfv
{
	template<typename T, unsigned int m_>
	class tensor1D
	{
		T tensor1D_[m_];
	public:

		tensor1D();
		~tensor1D();
		//Operator[] overload
		T& operator[](int i){ return tensor1D_[i]; }
		//Operator[] overload for constant
		const T& operator[](int i)const{ return tensor1D_[i]; }
		int size() const{ return m_; }

		const tensor1D<T, m_>& operator = (const tensor1D<T, m_>& m);
		tensor1D<T, m_> operator+(const tensor1D<T, m_>& m2) const;
		tensor1D<T, m_> operator-(const tensor1D<T, m_>& m2) const;
		tensor1D<T, m_> operator*(T a) const;
		tensor1D<T, m_> operator/(T a) const;

	};

	template<typename T, unsigned int m_>
	tensor1D<T, m_>::tensor1D()
	{
		memset(tensor1D_, 0, sizeof(T)*m_);
	}

	template<typename T, unsigned int m_>
	tensor1D<T, m_>::~tensor1D()
	{

	}

	template<typename T, unsigned int m_>
	tensor1D<T, m_> tensor1D<T, m_>::operator+(const tensor1D<T, m_>& m2)const
	{
		tensor1D<T, m_> result;
		for (int ii = 0; ii < m_; ++ii){
			result[ii] = (this->tensor1D_)[ii] + (m2.tensor1D_)[ii];
		}
		return result;
	}

	template<typename T, unsigned int m_>
	tensor1D<T, m_> tensor1D<T, m_>::operator-(const tensor1D<T, m_>& m2)const
	{
		tensor1D<T, m_> result;
		for (int ii = 0; ii < m_; ++ii){
			result[ii] = (this->tensor1D_)[ii] - (m2.tensor1D_)[ii];
		}
		return result;
	}

	template<typename T, unsigned int m_>
	tensor1D<T, m_> tensor1D<T, m_>::operator*(T a)const
	{
		tensor1D<T, m_> result;
		for (int ii = 0; ii < m_; ++ii){
			result[ii] = a * (this->tensor1D_)[ii];
		}
		return result;
	}


	template<typename T, unsigned int m_>
	tensor1D<T, m_> tensor1D<T, m_>::operator/(T a)const
	{
		tensor1D<T, m_> result;
		for (int ii = 0; ii < m_; ++ii){
			result[ii] = (this->tensor1D_)[ii] / a;
		}
		return result;
	}

	template<typename T, unsigned int m_>
	const tensor1D<T, m_>& tensor1D<T, m_>::operator = (const tensor1D<T, m_>& m){
		memcpy(tensor1D_, m.tensor1D_, sizeof(T)*m_);
		return *this;
	}

	template<typename T, unsigned int m_, unsigned int n_>
	class tensor2D
	{
	public:
		T tensor2D_[m_][n_];
		//tensor1D<T, m_> t1;
	public:
		tensor2D();
		~tensor2D();
		//tensor2D(const tensor2D<T, m_, n_>& MC);

		const std::pair<int, int> size()const { return std::pair<int, int>(m_, n_); }
		T* operator[](const int i);
		void print(std::ostream &out);
		const tensor2D<T, m_, n_>& operator = (const tensor2D<T, m_, n_>& m);
		tensor2D<T, m_, n_> operator+(const tensor2D<T, m_, n_>& m2) const;
		tensor2D<T, m_, n_> operator-(const tensor2D<T, m_, n_>& m2) const;
		tensor2D<T, m_, n_> operator*(T a) const;
		tensor2D<T, m_, n_> operator/(T a) const;
		//matrix dot matrix 
		template <unsigned int k_>
		tensor2D<T, m_, k_> operator*(const tensor2D<T, n_, k_>& m2) const;
	};

	//Constructor
	template<typename T, unsigned int m_, unsigned int n_>
	tensor2D<T, m_, n_>::tensor2D()
	{
		memset(*tensor2D_, 0, sizeof(T)*m_*n_);
	}

	//template<typename T, unsigned int m_, unsigned int n_>
	//tensor2D<T, m_, n_>::tensor2D(const tensor2D<T, m_, n_>& m)
	//{
	//	memcpy(*tensor2D_, *m.tensor2D_, sizeof(T)*n_*m_);
	//}

	template<typename T, unsigned int m_, unsigned int n_>
	tensor2D<T, m_, n_>::~tensor2D()
	{
	}

	template<typename T, unsigned int m_, unsigned int n_>
	T* tensor2D<T, m_, n_>::operator[](const int i)
	{
		return tensor2D_[i];//
	}

	template<typename T, unsigned int m_, unsigned int n_>
	const tensor2D<T, m_, n_>& tensor2D<T, m_, n_>::operator=(const tensor2D<T, m_, n_>& m)
	{
		memcpy(tensor2D_, m.tensor2D_, sizeof(T)*n_*m_);
		return *this;
	}


	template<typename T, unsigned int m_, unsigned n_>
	template <unsigned int k_>
	tensor2D<T, m_, k_> tensor2D<T, m_, n_>::operator*(const tensor2D<T, n_, k_>& m2)const
	{
		tensor2D<T, m_, k_> result;
		for (int ii = 0; ii < m_; ++ii){
			for (int jj = 0; jj < k_; ++jj){
				result[ii][jj] = 0.0;
				for (int kk = 0; kk < n_; ++kk){
					result[ii][jj] += (this->tensor2D_)[ii][kk] * (m2.tensor2D_)[kk][jj];
				}
			}
		}
		return result;
	}

	template<typename T, unsigned int m_, unsigned n_>
	tensor2D<T, m_, n_> tensor2D<T, m_, n_>::operator+(const tensor2D<T, m_, n_>& m2)const
	{
		tensor2D<T, m_, n_> result;
		for (int ii = 0; ii < m_; ++ii){
			for (int jj = 0; jj < n_; ++jj){
				result[ii][jj] = (this->tensor2D_)[ii][jj] + (m2.tensor2D_)[ii][jj];
			}
		}
		return result;
	}

	template<typename T, unsigned int m_, unsigned n_>
	tensor2D<T, m_, n_> tensor2D<T, m_, n_>::operator-(const tensor2D<T, m_, n_>& m2)const
	{
		tensor2D<T, m_, n_> result;
		for (int ii = 0; ii < m_; ++ii){
			for (int jj = 0; jj < n_; ++jj){
				result[ii][jj] = (this->tensor2D_)[ii][jj] - (m2.tensor2D_)[ii][jj];
			}
		}
		return result;
	}

	template<typename T, unsigned int m_, unsigned n_>
	tensor2D<T, m_, n_> tensor2D<T, m_, n_>::operator*(T a)const
	{
		tensor2D<T, m_, n_> result;
		for (int ii = 0; ii < m_; ++ii){
			for (int jj = 0; jj < n_; ++jj){
				result[ii][jj] = a * (this->tensor2D_)[ii][jj];
			}
		}
		return result;
	}

	template<typename T, unsigned int m_, unsigned n_>
	tensor2D<T, m_, n_> tensor2D<T, m_, n_>::operator/(T a)const
	{
		tensor2D<T, m_, n_> result;
		for (int ii = 0; ii < m_; ++ii){
			for (int jj = 0; jj < n_; ++jj){
				result[ii][jj] = (this->tensor2D_)[ii][jj] / a;
			}
		}
		return result;
	}

	template<typename T, unsigned int m_, unsigned int n_>
	void tensor2D<T, m_, n_>::print(std::ostream &out)
	{
		out << "{";
		for (int ii = 0; ii < m_; ++ii){
			out << "{";
			for (int jj = 0; jj < n_; ++jj){
				out << std::fixed << std::setw(20) << std::setprecision(10) << tensor2D_[ii][jj] << " ";
				if (jj != n_ - 1) out << ",";
			}
			out << "}";
			if (ii != m_ - 1) out << ",";
			out << std::endl;
		}
		out << "}";
	}

}

#endif
