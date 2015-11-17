#include <no_comment/DTK_DealIIEntityIntegrationRule.h>
#include <no_comment/DTK_DealIIHelpers.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

namespace DataTransferKit {

template <int dim>
void
getQuadratureRule(
    dealii::Quadrature<dim> const & dealii_quadrature,
    Teuchos::Array<Teuchos::Array<double> >& reference_points,
    Teuchos::Array<double>& weights )
{
    int const n = dealii_quadrature.size();
    reference_points.resize(n);
    weights.resize(n);
    for (int q = 0; q < n; ++q)
    {
        reference_points[q] = DealIIHelpers::getArray(dealii_quadrature.point(q));
        weights[q] = dealii_quadrature.weight(q);
    }
}



void
DealIIEntityIntegrationRule::
getIntegrationRule(
    const Entity& entity,
    const int order,
    Teuchos::Array<Teuchos::Array<double> >& reference_points,
    Teuchos::Array<double>& weights ) const
{
    int n = 0;
    while (2*n-1 < order)
        ++n;

    if (entity.topologicalDimension() == 3) {
        getQuadratureRule(
            dealii::QGauss<3>(n),
            reference_points,
            weights );
    } else if (entity.topologicalDimension() == 2) {
        getQuadratureRule(
            dealii::QGauss<2>(n),
            reference_points,
            weights );
    } else {
        throw std::runtime_error("EntityIntegrationRule not implemented");
    }
}

} // end namespace DataTransferKit
