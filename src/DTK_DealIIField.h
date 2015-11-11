#ifndef DTK_DEALIIFIELD_HPP
#define DTK_DEALIIFIELD_HPP

#include <no_comment/DTK_DealIITypes.h>
#include <DTK_Field.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <deal.II/dofs/dof_handler.h>

namespace DataTransferKit {

template <int dim,int spacedim>
class DealIIField : public Field
{
public:
    DealIIField(
        Teuchos::RCP<dealii::DoFHandler<dim,spacedim> const> const & dof_handler,
        Teuchos::RCP<DealIIVector> const & vector,
        unsigned int const & component );

    int dimension() const;

    Teuchos::ArrayView<const SupportId>
    getLocalSupportIds() const override;

    double readFieldData( const SupportId support_id,
        const int dimension ) const override;

    void writeFieldData( const SupportId support_id,
       const int dimension,
       const double data ) override;

    void finalizeAfterWrite() override;

  private:

    // Deal.II Vector.
    Teuchos::RCP<DealIIVector> d_dealii_vector;

    // Deal.II Degrees of Freedom Handler.
    Teuchos::RCP<dealii::DoFHandler<dim,spacedim> const> d_dealii_dof_handler;

    // Component in the system of equations.
    unsigned int d_component;

    // The support ids of the entities over which the field is constructed.
    Teuchos::Array<SupportId> d_support_ids;
};

} // end namespace DataTransferKit

#endif
