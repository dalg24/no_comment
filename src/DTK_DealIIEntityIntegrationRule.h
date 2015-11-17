#ifndef DTK_DEALIIENTITYINTEGRATIONRULE_HPP
#define DTK_DEALIIENTITYINTEGRATIONRULE_HPP

#include <DTK_EntityIntegrationRule.hpp>

namespace DataTransferKit {

class DealIIEntityIntegrationRule : public EntityIntegrationRule
{
public:

    void getIntegrationRule(
        const Entity& entity,
        const int order,
        Teuchos::Array<Teuchos::Array<double> >& reference_points,
        Teuchos::Array<double>& weights ) const override;

private:

};

} // end namespace DataTransferKit

#endif
