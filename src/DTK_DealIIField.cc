#include <no_comment/DTK_DealIIField.h>

namespace DataTransferKit {

DealIIField::
DealIIField(
    Teuchos::RCP<DealIIVector> const & vector )
    : d_dealii_vector(vector)
{
    dealii::IndexSet const local_indices =
        d_dealii_vector->locally_owned_elements();
    // TODO: ElementIterator (deal.II/base/index_set.h) does not follow std
    // input iterators specs
    // I am too lazy to fix it in deal.II so let's use a work around
//    d_support_ids.assign(local_indices.begin(), local_indices.end());
    d_support_ids.assign(local_indices.n_elements(), -1);
    auto it = d_support_ids.begin();
    for (auto dummy : local_indices)
        *(it++) = dummy;
    if (it != d_support_ids.end())
        throw std::runtime_error("something went wrong");
}

int
DealIIField::
dimension() const
{
    return 1;
}

Teuchos::ArrayView<const SupportId>
DealIIField::
getLocalSupportIds() const
{
    return d_support_ids();
}

double
DealIIField::
readFieldData(
    const SupportId support_id,
    const int dimension ) const
{
    if (dimension != 0)
        throw std::runtime_error("dimension should be zero...");
    return (*d_dealii_vector)[support_id];
}


void
DealIIField::
writeFieldData(
    const SupportId support_id,
    const int dimension,
    const double data)
{
    if (dimension != 0)
        throw std::runtime_error("dimension should be zero...");
    (*d_dealii_vector)[support_id] = data;
}


void
DealIIField::
finalizeAfterWrite()
{
    d_dealii_vector->update_ghost_values();
}

} // end namespace DataTransferKit
