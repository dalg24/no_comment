#ifndef DTK_DEALIINODALSHAPEFUNCTION_HPP
#define DTK_DEALIINODALSHAPEFUNCTION_HPP

#include <DTK_EntityShapeFunction.hpp>
#include <DTK_Types.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/fe/fe_values.h>
//#include <deal.II/base/quadrature.h>

namespace DataTransferKit
{

template <int dim,int spacedim>
class DealIIEntityShapeFunction : public EntityShapeFunction
{
public:

    DealIIEntityShapeFunction( 
      	const Teuchos::RCP<dealii::DoFHandler<dim,spacedim>>& dealii_dof_handler );

    /*!
     * \brief Given an entity, get the ids of its support locations.
     * \param entity Get the support locations for this entity.
     * \param support_ids Return the ids of the degrees of freedom in the parallel
     * vector space supporting the entities.
     */
    void entitySupportIds(
       	const Entity& entity,
       	Teuchos::Array<SupportId>& support_ids ) const;

    /*!
     * \brief Given an entity and a reference point, evaluate the shape
     * function of the entity at that point.
     * \param entity Evaluate the shape function of this entity.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param values Entity shape function evaluated at the reference
     * point. 
     */
    void evaluateValue( 
      	const Entity& entity,
      	const Teuchos::ArrayView<const double>& reference_point,
      	Teuchos::Array<double> & values ) const;

    /*!
     * \brief Given an entity and a reference point, evaluate the gradient of
     * the shape function of the entity at that point.
     * \param entity Evaluate the shape function of this entity.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param gradients Entity shape function gradients evaluated at the
     * reference point. Return these ordered with respect to those return by
     * getDOFIds() such that gradients[N][D] gives the gradient value of the
     * Nth DOF in the Dth spatial dimension.
     */
    void evaluateGradient( 
      	const Entity& entity,
        const Teuchos::ArrayView<const double>& reference_point,
      	Teuchos::Array<Teuchos::Array<double> >& gradients ) const;

private:

    // Deal.II Degrees of Freedom Handler.
    Teuchos::RCP<dealii::DoFHandler<dim,spacedim>> d_dealii_dof_handler;
//    // Deal.II Finite Element Values.
//    Teuchos::RCP<dealii::FEValues<dim,spacedim>>   d_dealii_fe_values;
//    // Deal.II Quadrature (points of interest where FE are evaluated)
//    Teuchos::RCP<dealii::Quadrature<dim>>          d_dealii_quadrature;

};

} // end namespace DataTransferKit

#endif
