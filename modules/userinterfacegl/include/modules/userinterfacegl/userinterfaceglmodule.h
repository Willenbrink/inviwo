/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2017-2018 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#ifndef IVW_USERINTERFACEGLMODULE_H
#define IVW_USERINTERFACEGLMODULE_H

#include <modules/userinterfacegl/userinterfaceglmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

#include <modules/userinterfacegl/glui/widgetsupplier.h>
#include <modules/userinterfacegl/glui/widgetfactory.h>
#include <modules/userinterfacegl/glui/widgetfactoryobject.h>
#include <modules/userinterfacegl/glui/element.h>

namespace inviwo {

class IVW_MODULE_USERINTERFACEGL_API UserInterfaceGLModule : public InviwoModule, public glui::WidgetSupplier {
public:
    UserInterfaceGLModule(InviwoApplication* app);
    virtual ~UserInterfaceGLModule();
    
    glui::WidgetFactory& getGLUIWidgetFactory();
    const glui::WidgetFactory& getGLUIWidgetFactory() const;

private:
    glui::WidgetFactory widgetFactory_;
};

}  // namespace inviwo

#endif  // IVW_USERINTERFACEGLMODULE_H
