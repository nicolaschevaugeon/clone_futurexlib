#ifndef small_any_h
#define small_any_h

#include <type_traits> 
#include <stdio.h>
#include <string.h>


namespace xinterface{

  namespace xmeshinterface{


    /*
      1 char                    =  1 byte  = 8 bits
      1 short                   =  2 bytes
      1 int                     =  4 bytes (1 int = 2 short)
      1 long                    =  8 bytes (1 long = 2 int)  
      1 long long               =  8 bytes
      1 void*                   =  4 bytes (1 void* = 1 int)
      1 pair<int,int>           =  8 bytes ( = 2 int)
      1 pair<int,pair<int,int>> = 12 bytes ( = 3 int)
    */

#define MAXSIZE  4
  class small_any;
    template<class T>
    inline const typename std::remove_reference<T>::type &any_cast( const small_any & any);
    template<class T>
    inline typename std::remove_reference<T>::type &any_cast( small_any & any);

    class small_any{
    public:
    private:
      using storage_type =  std::aligned_storage<MAXSIZE*sizeof(void *), MAXSIZE*sizeof(void *)>::type;
      // public:
      storage_type buffer;
    public:

      template <class T>
      friend inline const typename std::remove_reference<T>::type &any_cast( const small_any & any);
      template<class T>
      friend inline typename std::remove_reference<T>::type &any_cast( small_any & any);

      template<class T>
	small_any( const T & in)
	{
	  static_assert ( (sizeof(T) <= sizeof(storage_type) ), "The object can not be contained in small_any");

          any_cast<T>(*this) = in;
          //(typename std::remove_reference<T>::type &)(this->buffer) = in;
        }

      small_any( ) =default;


      small_any(const small_any &other) =default;
      small_any& operator=( const small_any & ) = default;
      small_any& operator=( small_any && ) = default;
      small_any(small_any &&other):buffer(std::move(other.buffer))    {    }

      template<class T>
	small_any& operator=( const T & in)
	{
	  static_assert ( (sizeof(T) <= sizeof(storage_type) ), "The object can not be contained in small_any");
	  memcpy(&buffer, &in, sizeof(T));
	  return *this;
        }

      inline friend bool operator==(const small_any& left, const small_any& right){ 
	return (memcmp( &left.buffer, &right.buffer , sizeof(storage_type)  ) ==0 );
      }
      /*
	inline friend bool operator<(const small_any& left, const small_any& right){ 
	return (memcmp( &left.buffer, &right.buffer , sizeof(storage_type)  ) <0 );
	}

	inline friend bool operator>(const small_any& left, const small_any& right){ 
	return (memcmp( &left.buffer, &right.buffer , sizeof(storage_type)  ) >0 );
	}
      */
    };

    //C++11
    template<class T>
      inline typename std::remove_reference<T>::type &any_cast( small_any & any)
      {
        return (typename std::remove_reference<T>::type &)(any.buffer);
      }
    template<class T>
      inline const typename std::remove_reference<T>::type &any_cast( const small_any & any)
      {
        return (typename std::remove_reference<T>::type &)(any.buffer);
      }
      //C++14
      /*
        template<class T>
        inline const std::remove_reference_t<T> &any_cast( const small_any & any){
        return (std::remove_reference_t<T> &)(any.buffer);
        }
        template<class T>
        inline std::remove_reference_t<T> &any_cast( small_any & any){
        return (std::remove_reference_t<T> &)(any.buffer);
        }
      */

  } // namepsace xmeshinterface
} // namespace xinterface

#endif
