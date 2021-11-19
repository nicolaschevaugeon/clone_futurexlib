#ifndef _X_FIELD_H
#error Do NOT include xField_imp.h alone
#endif




// Start definition of xField Template members
template<typename VT>
template<class T>
xField<VT>::xField( xValueManagerDist<VT>* vm, const T& f) : value_manager(vm)
{
    _insert(spacePtr(new T(f)));
    newstorage();
}

template<typename VT>
template <class T1, class T2>
xField<VT>::xField(xValueManagerDist<VT>* vm, const T1& f1, const T2& f2) : value_manager(vm)
{
    _insert(spacePtr(new T1(f1)));
    _insert(spacePtr(new T2(f2)));
    newstorage();
}


template<typename VT>
template <class T1, class T2, class T3>
xField<VT>::xField(xValueManagerDist<VT>* vm, const T1& f1, const T2& f2, const T3& f3) : value_manager(vm)
{
    _insert(spacePtr(new T1(f1)));
    _insert(spacePtr(new T2(f2)));
    _insert(spacePtr(new T3(f3)));
    newstorage();
}

template<typename VT>
template <class S>
void xField<VT>::ResetStoragePolicy()
{
    delete storage;
    storage=new S(*value_manager,spaces);
}

template<typename VT>
template <class T>
void  xField<VT>::insert(const T& f)
{
    _insert(spacePtr(new T(f)));
    clearstorage();
}

template<typename VT>
template<class T>
T&  xField<VT>::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& v) const
{
    AOMD::mEntity* e = geo_appro->getEntity();
    getVal(beginFcts(e), endFcts(e), beginValues(e), endValues(e), sizeFcts(e), geo_appro, geo_integ, v);
    return v;
}

template<typename VT>
template<class T>
T&  xField<VT>::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& v) const
{
    AOMD::mEntity* e = geo_appro->getEntity();
    getGrad(beginFcts(e), endFcts(e), beginValues(e),endValues(e), sizeFcts(e),geo_appro, geo_integ, v);
    return v;
}

template<typename VT>
template<class T>
T&   xField<VT>::getGradLocal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, T& v) const
{
    AOMD::mEntity* e = geo_appro->getEntity();
    getGradLocal(beginFcts(e), endFcts(e), beginValues(e), endValues(e), geo_appro, geo_integ, v);
    return v;
}

template<typename VT>
template <class iterFct, class UnaryOperator>
void xField<VT>::getFF(iterFct it, iterFct ite,
                       std::vector<typename UnaryOperator::result_type>& ff,
                       const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                       UnaryOperator funct, bool resetStorage)
{
    typename UnaryOperator::argument_type v;
    for (; it != ite; ++it)
    {
        if(resetStorage){(*it)->resetStorage();}
        (*it)->getVal(geo_appro, geo_integ, v);
        //          cout<<v<<endl;
        ff.push_back(funct(v));
    }
}

template<typename VT>
template <class iterFct, class UnaryOperator>
void xField<VT>::getGradFF (iterFct it, iterFct ite,
                            std::vector<typename  UnaryOperator::result_type >& grads,
                            const xGeomElem* geo_appro,const xGeomElem* geo_integ,
                            UnaryOperator funct, bool resetStorage)
{
    typename UnaryOperator::argument_type v;
    for (; it != ite; ++it)
    {
        if(resetStorage){(*it)->resetStorage();}
        (*it)->getGrad(geo_appro, geo_integ, v);
        grads.push_back(funct(v));
    }
}

template<typename VT>
template <class iterFct, class UnaryOperator, class CONTAINER>
void xField<VT>::getGradFFEigen (iterFct it, iterFct ite,
                                 CONTAINER& grads,
                                 const xGeomElem* geo_appro,const xGeomElem* geo_integ,
                                 UnaryOperator funct, bool resetStorage)
{
    typename UnaryOperator::argument_type v;
    for (int i=0; it != ite; ++it, ++i)
    {
        if(resetStorage){(*it)->resetStorage();}
        (*it)->getGrad(geo_appro, geo_integ, v);
        grads.col(i)=funct(v);
    }
}

template<typename VT>
template <class iterFct, class UnaryOperator>
void xField<VT>::getGradLocalFF (iterFct it, iterFct ite,
                                 std::vector<typename  UnaryOperator::result_type >& grads,
                                 const xGeomElem* geo_appro,const xGeomElem* geo_integ,
                                 UnaryOperator funct)
{
    typename UnaryOperator::argument_type v;
    for (; it != ite; ++it)
    {
        (*it)->getGradLocal(geo_appro, geo_integ, v);
        grads.push_back(funct(v));
    }
}

template<typename VT>
template <class iterFct>
void xField<VT>::resetStorage (iterFct it, iterFct ite,
                               const xGeomElem* geo_appro, const xGeomElem* geo_integ)
{
    for (; it != ite; ++it) (*it)->resetStorage();
}


template<typename VT>
template <class iterFct, class iterVal, class T>
void xField<VT>::getVal (iterFct it, iterFct ite,
                                iterVal first, iterVal last,
                                int size_Fct_Val,
                                const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                T& v)
{
    std::vector<T> ff;
    ff.reserve(size_Fct_Val);
    getFF(it, ite, ff, geo_appro, geo_integ, xtool::xIdentity<T>());//-----
    v = inner_product(ff.begin(), ff.end(), first, T(), std::plus<T>(), ptr_product<T>());
    //   cout<<v<<endl;
    return;
}


template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getVal (iterFct it, iterFct ite,
                                iterVal first, iterVal last,
                                int size_Fct_Val,
                                const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                std::complex<double>& v)
{
    std::vector<double > ff;
    ff.reserve(size_Fct_Val);
    getFF(it, ite, ff, geo_appro, geo_integ, xtool::xIdentity<double>());//-----
    v = inner_product(ff.begin(), ff.end(), first, std::complex<double>(), std::plus<std::complex<double> >(), ptr_product<std::complex<double> >());
    //   cout<<v<<endl;
    return;
}



template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getVal (iterFct it, iterFct ite,
                                iterVal first, iterVal last,
                                int size_Fct_Val,
                                const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                xtensor::xVectorDoubleComplex& v)
{
    std::vector<xtensor::xVector<double> > ff;
    ff.reserve(size_Fct_Val);
    getFF(it, ite, ff, geo_appro, geo_integ, xtool::xIdentity<xtensor::xVector<double> >());//-----
    v = inner_product(ff.begin(), ff.end(), first, xtensor::xVectorDoubleComplex(), std::plus<xtensor::xVectorDoubleComplex >(), ptr_product<xtensor::xVectorDoubleComplex >());
    //   cout<<v<<endl;
    return;
}


template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getGrad(iterFct it, iterFct ite,
                                iterVal first, iterVal last,
                                int size_Fct_Val,
                                const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                xtensor::xVector<>& v)
{
    std::vector<xtensor::xVector<> > grads;
    grads.reserve(size_Fct_Val);
    getGradFF(it, ite, grads, geo_appro, geo_integ, xtool::xIdentity<xtensor::xVector<> >() );
    v = inner_product(grads.begin(), grads.end(), first, xtensor::xVector<>(), std::plus<xtensor::xVector<> >(), ptr_product<xtensor::xVector<> >());
    return;
}


template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getGrad(iterFct it, iterFct ite,
                                iterVal first, iterVal last,
                                int size_Fct_Val,
                                const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                xtensor::xVectorDoubleComplex& v)
{
    // In order to simplify inner_product implementation
    // we cast xVector<double> gradients to xVector<complex<double> >
    // Note that we could also proceed by hand, just like for xTensor2<double>
    std::vector<xtensor::xVectorDoubleComplex > grads;
    grads.reserve(size_Fct_Val);
    getGradFF(it, ite, grads, geo_appro, geo_integ, xtensor::xCastToVectorComplex<double>() );// Note the cast
    v = inner_product(grads.begin(), grads.end(), first, xtensor::xVectorDoubleComplex(), std::plus<xtensor::xVectorDoubleComplex >(), ptr_product<xtensor::xVectorDoubleComplex >());

    return;
}


template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getGrad
(iterFct it, iterFct ite,
 iterVal first, iterVal last,
 int size_Fct_Val,
 const xGeomElem* geo_appro, const xGeomElem* geo_integ,
 xtensor::xTensor2<>& v)
{
    std::vector<xtensor::xTensor2<> > grads;
    grads.reserve(size_Fct_Val);
    getGradFF(it, ite, grads, geo_appro, geo_integ, xtool::xIdentity<xtensor::xTensor2<> >() );
    std::vector<xtensor::xTensor2<> >::iterator itt = grads.begin();
    for (; itt != grads.end(); ++itt, first++)
    {
        xtensor::xTensor2<>& g = *itt;
        axpy_tensor2((*first)->getVal(), g, v);
    }
    return;
}



template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getGrad(iterFct it, iterFct ite,
                                iterVal first, iterVal last,
                                int size_Fct_Val,
                                const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                xtensor::xTensor2DoubleComplex& v)
{
    std::vector<xtensor::xTensor2DoubleComplex > grads;
    grads.reserve(size_Fct_Val);
    getGradFF(it, ite, grads, geo_appro, geo_integ, xtensor::xCastToTensor2Complex<double>() );
    auto itt = grads.begin();
    for (; itt != grads.end(); ++itt, first++)
    {
        xtensor::xTensor2DoubleComplex& g = *itt;
        axpy_tensor2((*first)->getVal(), g, v);
    }
    return;
}



template<typename VT>
template <class iterFct, class iterVal>
void xField<VT>::getGradLocal (iterFct it, iterFct ite,
                               iterVal first, iterVal last,
                               const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                               xtensor::xTensor2<>& v)
{
    for(;first != last; first++, it++)
    {
        VT val = (*first)->getVal();
        (*it)->getGradLocalAxpy(geo_appro, geo_integ, v, val);
    }
    return;
}

// End implementation of xField Template members


//-------------------------------------------------------------------------
template<typename VT>
xField<VT>::~xField(){
    delete storage;
}


template<typename VT>
xField<VT>::xField(xValueManagerDist<VT>* vm) : value_manager(vm)
{
    newstorage();
}

template<typename VT>
xField<VT>::xField(xValueManagerDist<VT>* vm, const spacePtr& f) : value_manager(vm)
{
    _insert(f);
    newstorage();
}

template<typename VT>
void xField<VT>::newstorage(){
    storage = new xFieldStorageElement<VT>(*value_manager, spaces);
}

template<typename VT>
xFieldStorage<VT> *  xField<VT>::GetStoragePolicy()
{
    return storage;
}

template<typename VT>
void xField<VT>::ResetStoragePolicy()
{
    delete storage;
    newstorage();
}

template<typename VT>
void xField<VT>::clearstorage(){
    storage->clear();
}

template<typename VT>
void xField<VT>::insert(const spacePtr& f)
{
    _insert(f);
    clearstorage();
}

template<typename VT>
int xField<VT>::sizeFcts(AOMD::mEntity* e) const
{
    return storage->sizeFcts(e);
}

template<typename VT>
typename std::vector<typename xField<VT>::shapeFctPtr>::iterator xField<VT>::beginFcts(AOMD::mEntity* e) const
{
    return storage->beginFcts(e);
}

template<typename VT>
typename std::vector<typename xField<VT>::shapeFctPtr>::iterator xField<VT>::endFcts(AOMD::mEntity* e) const
{
    return storage->endFcts(e);
}

template<typename VT>
auto xField<VT>::beginValues(AOMD::mEntity* e)  const -> typename value_container_t::iterator
{
    return storage->beginValues(e);
}

template<typename VT>
auto xField<VT>::endValues(AOMD::mEntity* e) const -> typename value_container_t::iterator
{
    return storage->endValues(e);
}

template<typename VT>
void  xField<VT>::clear()    {
    spaces.clear();
    clearstorage();
}


template<typename VT>
void xField<VT>::_insert(spacePtr d)
{
    //  if (spaces.empty()) TensoDim = d->getTensoDim();
    //  else assert(TensoDim == d->getTensoDim());
    spaces.push_back(d);
}

// void xField::update(mEntity* e) const
// {
//   const bool debug = xdebug_flag;
//   if (e != fem_curr.getEntity() )
//     {
//       fem_curr.setKeysAndFcts(e, begin(), end());
//       vals_curr_ptr.clear();
//       value_manager->getValPtr(fem_curr.beginKey(), fem_curr.endKey(), vals_curr_ptr);
//       if (debug) std::cout << " in update " << std::endl;
//     }
// }

template<typename VT>
void xField<VT>::getVals(AOMD::mEntity* e, std::vector<VT>& vals)
{
    //update(e);
    //iter_ptr it = vals_curr_ptr.begin(), itEnd = vals_curr_ptr.end();
    //for ( ; it != itEnd; ++it) vals.push_back((*it)->getVal());

    iter_ptr it = storage->beginValues(e), itEnd = storage->endValues(e);
    for ( ; it != itEnd; ++it) vals.push_back((*it)->getVal());

}

template<typename VT>
void xField<VT>::setVal(AOMD::mEntity* e, const VT& v)
{
    //update(e);
    for (auto it =
         storage->beginValues(e); it != storage->endValues(e); ++it)
    {
        (*it)->setVal(v);
    }

    return;
}

template<typename VT>
xFillFieldFromZeroForm<VT>::xFillFieldFromZeroForm(xField<VT>& f, xFormZero<VT>& fo)
    :  xIntegrateFormCommand(&fo), field(f), form(&fo), total(0.)   {}


template<typename VT>
void  xFillFieldFromZeroForm<VT>::openApproxElem(xGeomElem* g_appro)
{
    geom_appro = g_appro;
    form->init(g_appro);
}

template<typename VT>
void  xFillFieldFromZeroForm<VT>::closeApproxElem(xGeomElem* g_appro)
{
    total += form->getScalar().getVal();
    AOMD::mEntity* e_integ = geom_integ->getEntity();
    field.setVal(e_integ, form->getScalar().getVal());
}


//-------------------------------------------------------------------------
